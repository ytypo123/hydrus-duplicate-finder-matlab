%% Hydrus Duplicate Finder & Potential Duplicate Setter
%
% ------------------------------------------------------
% QUICK START (for Hydrus users)
% ------------------------------------------------------
% 1. In Hydrus Client, enable "Client API"
% 2. Get an API key:
% 3. Create a configuration file named "hydrus_config.json" in the same
%    folder as this script. Example:
%       {
%         "api_key": "your-hydrus-api-key-here",
%         "api_url": "http://127.0.0.1",
%         "api_port": 42001,
%       }
%
% 4. Run in MATLAB. The script will:
%       - Load settings from hydrus_config.json
%       - Fetch file hashes + metadata from Hydrus
%       - Build candidate pairs (duration, AR)
%       - Optional: perceptual-hash prefilter (aHash)
%       - Compare thumbnails (SSIM/MS-SSIM)
%       - Post "potential duplicate" relationships
%       - Print top tags among duplicate files and save a full CSV
%
% Notes:
%   • Only posts "potential duplicate" (relationship = 0) to Hydrus.
%   • Already-marked pairs may return "400 Bad Request" — counted as skipped.
%   • MATLAB R2021a+ recommended for MS-SSIM (`multissim` function).
%
% Requirements:
%   - MATLAB + Image Processing Toolbox
%   - Hydrus Client API running and accessible
% ------------------------------------------------------

%% -------------------- CONFIGURATION --------------------
clc
clearvars
close all
tAll = tic;

% Load configuration from JSON file
configFile = 'hydrus_config.json';
if ~isfile(configFile)
    error('Configuration file "%s" not found.', configFile);
end

fid = fopen(configFile, 'r');
raw = fread(fid, inf);
fclose(fid);

config = jsondecode(char(raw'));

% Extract API details
api_key = config.api_key;
api_URL = sprintf('%s:%d', config.api_url, config.api_port);

% Tags to search for in Hydrus
tags = {'system:filetype is video', 'system:filesize > 700MB'};

% Verbose logging (true = detailed debug prints inside helpers)
VERBOSE = false;

% -------- Filter toggles --------
USE_AR_FILTER       = true;   % If false, skip aspect ratio filter
USE_AHASH_PREFILTER = true;   % If false, skip aHash prefilter (faster, but more SSIM work)

% Candidate generation tolerances
tolSeconds = 1.0;    
% Duration tolerance (seconds) for initial candidate filtering.
% Example: 2.0 = ±2 s, 10.0 = ±10 s. Larger => more pairs => slower SSIM step.
% Note: Real duplicates can differ by tens of seconds (e.g., padding/credits).

% Aspect ratio tolerance (fractional). 0.01 = ±1% difference allowed.
arTol = 0.01;
% Example for 16:9 AR ≈ 1.7778:
%     0.005 → [1.769, 1.787] (very strict)
%     0.010 → [1.760, 1.796] (common)
%     0.020 → [1.742, 1.813] (looser)
%     0.050 → [1.689, 1.867] (very loose)

% SSIM / preprocessing parameters
P = struct();
P.fixedResize   = [256 256];   % Smaller = faster, larger = more accurate
P.toGray        = true;        % Ignore color; compare structure/texture
P.preBlurSigma  = 0.6;         % Small blur (0.4–1.0) to reduce artifacts
P.histEq        = false;       % Equalize histograms (useful if lighting varies)
P.useMSSSIM     = true;        % Use MS-SSIM if available (R2021a+)
P.ssimK         = [0.01 0.03]; % SSIM constants (leave default)
P.ssimExponents = [1 1 1];     % MS-SSIM exponents (leave default)

% Similarity threshold for declaring duplicates (0..1).
similarity_threshold = 0.50;   % Lower finds more (risk FP), higher is stricter.

% Posting options for setFileRelationships
postOpts = struct();
postOpts.batch_size       = 200;  % Hydrus batch size
postOpts.max_retries      = 3;    % Retry attempts for network errors
postOpts.retry_pause_s    = 1.0;  % Wait between retries
postOpts.timeout_s        = 60;   % Web request timeout
postOpts.do_default_merge = true; % Let Hydrus auto-merge tags etc.

% -------- Tag console report controls --------
TopKConsoleTags  = 10;   % Show at most N tags in console
MinConsoleCount = 2;   % show only tags that appear in >= 2 duplicate files

%% -------------------- STEP 0: Sanity check API config --------------------
% Fail fast with friendly message if API details are missing/malformed.
if ~(ischar(api_key) || isstring(api_key)) || strlength(string(api_key)) < 10
    error(['Hydrus API key is missing or too short.\n' ...
           'Go to Hydrus: Services -> Manage Services -> Client API -> Permissions -> API Key,\n' ...
           'copy a key, and paste it into api_key at the top of this script.']);
end
if ~(ischar(api_URL) || isstring(api_URL)) || ~startsWith(string(api_URL), ["http://","https://"])
    error('api_URL must start with http:// or https:// (e.g., http://127.0.0.1:42001)');
end

%% -------------------- STEP 1: Validate Hydrus connection --------------------
validate_connection(api_key, api_URL, VERBOSE);

%% -------------------- STEP 2: Fetch file hashes --------------------
t = tic;
if VERBOSE
    fprintf('Getting these tags: %s\n', strjoin(tags, ', '));
end
file_hash = get_hashes(tags, api_key, api_URL, VERBOSE);
if isstring(file_hash), file_hash = cellstr(file_hash); end
n = numel(file_hash);
fprintf('Got %d hashes in (%.2f s)\n', n, toc(t));
if n == 0
    fprintf('No files found for given tags. Exiting.\n');
    return
end

%% -------------------- STEP 3: Fetch metadata --------------------
t = tic;
[file_size, file_duration, file_num_frames, file_width, file_height, file_aspect_ratio] = ...
    get_metadata(file_hash, api_key, api_URL, VERBOSE);
fprintf('Fetched metadata for %d files in (%.2f s)\n', n, toc(t));

%% -------------------- STEP 4: Candidate generation --------------------
% 4.1 Duration filter
t = tic;
file_duration_pairs = findPairsByDuration(file_duration, tolSeconds);
fprintf('Duration candidates: %d pairs (%.2f s)\n', size(file_duration_pairs,1), toc(t));

% 4.2 Aspect ratio filter (optional)
if USE_AR_FILTER
    t = tic;
    file_aspect_ratio_pairs = findPairsByAR_tol(file_aspect_ratio, arTol, file_duration_pairs, true);
    fprintf('AR-filtered candidates: %d pairs (%.2f s)\n', size(file_aspect_ratio_pairs,1), toc(t));
else
    file_aspect_ratio_pairs = file_duration_pairs;
    fprintf('AR filter disabled: passing through %d duration-based pairs.\n', size(file_aspect_ratio_pairs,1));
end

if isempty(file_aspect_ratio_pairs)
    fprintf('No candidates after filtering. Exiting.\n');
    fprintf('Total time: %.1f s\n', toc(tAll));
    return
end

%% -------------------- STEP 5.0: Fast perceptual-hash prefilter (optional) --------------------
% Keeps only pairs whose perceptual hashes are close (small Hamming distance).
% Much cheaper than SSIM and can dramatically reduce SSIM workload.
if USE_AHASH_PREFILTER
    t = tic;
    aHashOpts = struct();
    aHashOpts.resize     = [64 64];   % tiny resize before 8x8 mean (fast)
    aHashOpts.hammingMax = 10;        % <= 10/64 bits (~15%); tune 6–12

    origCount = size(file_aspect_ratio_pairs,1);
    [pfPairs, pfStats] = prefilterPairsWithAHash( ...
        file_aspect_ratio_pairs, file_hash, api_key, api_URL, aHashOpts, VERBOSE);

    keptCount = size(pfPairs,1);
    pctKept   = 100 * keptCount / max(1, origCount);
    fprintf('aHash prefilter: kept %d / %d candidate pairs (%.1f%%) in (%.2f s)\n', ...
        keptCount, origCount, pctKept, toc(t));

    if pfStats.missing > 0
        fprintf('  (Skipped %d pairs due to missing/unreadable thumbnails)\n', pfStats.missing);
    end
    if pctKept < 5
        fprintf('  NOTE: Very aggressive filtering — if you miss matches, raise hammingMax a bit.\n');
    elseif pctKept > 50
        fprintf('  NOTE: Weak filtering — lower hammingMax to speed up SSIM.\n');
    end

    if isempty(pfPairs)
        fprintf('No pairs survived the perceptual-hash prefilter. Exiting.\n');
        fprintf('Total time: %.1f s\n', toc(tAll));
        return
    end

    % Use reduced set for SSIM
    file_aspect_ratio_pairs = pfPairs;
else
    fprintf('aHash prefilter disabled: passing through %d pairs to SSIM.\n', size(file_aspect_ratio_pairs,1));
end

%% -------------------- STEP 5: Thumbnail comparison (SSIM / MS-SSIM) --------------------
t = tic;
duplicatesMatrix = findAndCompareThumbnails( ...
    file_aspect_ratio_pairs, file_hash, api_key, api_URL, VERBOSE, similarity_threshold, P);
fprintf('Thumbnail compare found %d duplicate pairs in (%.2f s)\n', size(duplicatesMatrix,1), toc(t));

if isempty(duplicatesMatrix)
    fprintf('No duplicates above threshold %.3f. Nothing to post.\n', similarity_threshold);
    fprintf('Total time: %.1f s (%.1f files/s)\n', toc(tAll), n / toc(tAll));
    return
end

%% -------------------- STEP 6: Post potential duplicates to Hydrus --------------------
pairsToPost = size(duplicatesMatrix, 1);

fprintf('\n%s\n', repmat('=',1,72));
fprintf('FINAL RESULT\n');
fprintf('%s\n', repmat('-',1,72));
fprintf('Pairs to post: %d\n', pairsToPost);

t = tic;
rep = setFileRelationships(api_key, api_URL, file_hash, duplicatesMatrix, postOpts);
elapsed = toc(t);

fprintf('Relationships -> Posted: %d | Skipped: %d | Failed: %d | Attempted: %d\n', ...
        rep.posted, rep.skipped, rep.failed, rep.attempted);
fprintf('Elapsed time: %.2f s\n', elapsed);
fprintf('%s\n', repmat('-',1,72));
fprintf(['Note: "Posted" means the potential duplicate was sent to Hydrus.\n' ...
         '      If the pairs are later marks it "Not a Duplicate" in Hydrus, this value can still be > 0.\n']);
fprintf('%s\n\n', repmat('=',1,72));


%% -------------------- STEP 6.5: Tag report among duplicate files --------------------
% Console: show only the TOP 10 tags that appear in at least MinConsoleCount
% duplicate files. The full list (with counts & percentages) still goes to CSV.

% If you haven't added this to CONFIGURATION, we'll default to 2.
if ~exist('MinConsoleCount','var') || isempty(MinConsoleCount)
    MinConsoleCount = 2;  % e.g., show tags present in >= 2 duplicate files
end
if ~exist('TopKConsoleTags','var') || isempty(TopKConsoleTags)
    TopKConsoleTags = 10;
end

try
    dupIdx    = unique(duplicatesMatrix(:));
    dupHashes = file_hash(dupIdx);

    tagLists = get_tags_for_hashes(dupHashes, api_key, api_URL, VERBOSE);

    % (Optional) peek at one file's tags when debugging
    if VERBOSE
        nonEmpty = find(~cellfun(@isempty, tagLists), 1, 'first');
        if ~isempty(nonEmpty)
            fprintf('Example tags from one duplicate file:\n');
            disp(tagLists{nonEmpty}(1:min(end,20)));
        else
            fprintf('No tags parsed for duplicate files (all empty).\n');
        end
    end

    % Count all tags across duplicate files
    [allTags, allCounts] = countAllTags(tagLists);
    dupFileCount = numel(dupHashes);
    pct = (allCounts(:) ./ max(1, dupFileCount)) * 100;

    % Console: top K with count >= MinConsoleCount
    mask = allCounts >= MinConsoleCount;
    topTags  = allTags(mask);
    topCount = allCounts(mask);
    topPct   = pct(mask);

    if ~isempty(topTags)
        K = min(TopKConsoleTags, numel(topTags));
        fprintf('\n--- Top %d Tags Among Duplicate Files (Count ≥ %d) ---\n', K, MinConsoleCount);
        for k = 1:K
            fprintf('  %5d  (%4.1f%%%%)  %s\n', topCount(k), topPct(k), topTags{k});
        end
        fprintf('------------------------------------------------------\n\n');
    else
        fprintf('\nNo tags met the Count ≥ %d threshold for console display. See CSV for full list.\n\n', MinConsoleCount);
    end

    % CSV: full list with counts + percentages (no filtering)
    T = table(allCounts(:), pct, string(allTags(:)), ...
              'VariableNames', {'Count', 'PercentOfDupFiles', 'Tag'});
    ts = datestr(now, 'yyyy-mm-dd_HHMMSS');
    csvName = sprintf('hydrus_duplicate_tags_%s.csv', ts);
    writetable(T, csvName);
    fprintf('Saved full tag list to %s\n', csvName);

catch e
    fprintf('Full-tag report skipped due to error: %s\n', e.message);
end


%% -------------------- STEP 7: Final stats --------------------
elapsed = toc(tAll);
fprintf('Finished processing %d files in %.1f seconds (%.1f files/s)\n', n, elapsed, n/elapsed);

%% -------------------- HELPER FUNCTIONS --------------------
% All helper functions are embedded below so this script is 100% portable.

% --- validate_connection ---
function validate_connection(api_key, api_URL, VERBOSE)
    if VERBOSE
        disp('Testing if MATLAB can access Hydrus.');
    end
    try
        url = sprintf('%s/api_version', api_URL);
        options = weboptions('RequestMethod', 'GET');
        response = webread(url, options);
        if VERBOSE
            fprintf('Hydrus API Version number: %d\n', response.version);
            disp('MATLAB seems to be able to access Hydrus. Probably.');
        end
    catch err
        error('Failed to connect to Hydrus API: %s', err.message);
    end
end

% --- get_hashes ---
function hashes = get_hashes(tags, api_key, api_URL, VERBOSE)
    if VERBOSE
        disp({'Getting these tags', tags{:}});
    end
    searchEndpoint = '/get_files/search_files';
    requestURL = sprintf('%s%s?tags=%s&return_hashes=true&return_file_ids=false&Hydrus-Client-API-Access-Key=%s', ...
        api_URL, searchEndpoint, urlencode(jsonencode(tags)), api_key);
    options = weboptions('RequestMethod', 'GET');
    response = webread(requestURL, options);
    hashes = response.hashes;
    if VERBOSE
        fprintf('Got %d hashes\n', numel(hashes));
    end
end

% --- get_metadata ---
function [sizes, durations, num_frames, file_widths, file_heights, file_aspect_ratios] = ...
         get_metadata(hashes, api_key, api_URL, VERBOSE)
    searchEndpoint = '/get_files/file_metadata';
    batch_size = 100;
    n = numel(hashes);
    batches = fix(n / batch_size);
    remainder = mod(n, batch_size);
    sizes = zeros(1, n);
    durations = zeros(1, n);
    num_frames = zeros(1, n);
    file_widths = zeros(1, n);
    file_heights = zeros(1, n);
    options = weboptions('RequestMethod', 'GET');
    baseRequestURL = sprintf('%s%s?only_return_basic_information=true&Hydrus-Client-API-Access-Key=%s', ...
        api_URL, searchEndpoint, api_key);
    for i = 1:batches
        requestURL = sprintf('%s&hashes=%s', baseRequestURL, urlencode(jsonencode(hashes((i-1)*batch_size+1:i*batch_size))));
        responses = webread(requestURL, options);
        for j = 1:batch_size
            idx = (i-1)*batch_size + j;
            sizes(idx) = responses.metadata(j).size;
            durations(idx) = responses.metadata(j).duration;
            num_frames(idx) = responses.metadata(j).num_frames;
            file_widths(idx) = responses.metadata(j).width;
            file_heights(idx) = responses.metadata(j).height;
        end
        if VERBOSE
            fprintf('Processed batch %d/%d (%d items)\n', i, batches, batch_size);
        end
    end
    if remainder > 0
        requestURL = sprintf('%s&hashes=%s', baseRequestURL, urlencode(jsonencode(hashes(batches*batch_size+1:end))));
        responses = webread(requestURL, options);
        for j = 1:remainder
            idx = batches*batch_size + j;
            sizes(idx) = responses.metadata(j).size;
            durations(idx) = responses.metadata(j).duration;
            num_frames(idx) = responses.metadata(j).num_frames;
            file_widths(idx) = responses.metadata(j).width;
            file_heights(idx) = responses.metadata(j).height;
        end
    end
    file_aspect_ratios = file_widths ./ file_heights;
    if VERBOSE
        disp('Got all metadata');
    end
end

% --- findPairsByDuration ---
function pairs = findPairsByDuration(durations, tolSeconds)
    % Returns index pairs where durations differ by <= tolSeconds
    durations = durations / 1000; 
    n = numel(durations);
    pairs = [];
    for i = 1:n-1
        for j = i+1:n
            if abs(durations(i) - durations(j)) <= tolSeconds
                pairs = [pairs; i j]; %#ok<AGROW>
            end
        end
    end
end

% --- findPairsByAR_tol ---
function pairs_out = findPairsByAR_tol(aspect_ratios, tol, candidatePairs, relative)
    % Filters candidate pairs by aspect ratio difference
    if relative
        matchFn = @(ar1, ar2) abs(ar1 - ar2) <= tol * max(ar1, ar2);
    else
        matchFn = @(ar1, ar2) abs(ar1 - ar2) <= tol;
    end
    pairs_out = [];
    for k = 1:size(candidatePairs,1)
        i = candidatePairs(k,1);
        j = candidatePairs(k,2);
        if matchFn(aspect_ratios(i), aspect_ratios(j))
            pairs_out = [pairs_out; i j]; %#ok<AGROW>
        end
    end
end

% --- findAndCompareThumbnails ---
function duplicatesMatrix = findAndCompareThumbnails(potentialDuplicatePairs, fileHashes, apiKey, apiURL, verbose, similarityThreshold, P)
    % Compares candidate pairs of files by thumbnail similarity
    duplicatesMatrix = [];
    searchEndpoint = '/get_files/thumbnail';
    options = weboptions('RequestMethod', 'GET', 'ContentType', 'image');
    for i = 1:size(potentialDuplicatePairs, 1)
        requestURL1 = sprintf('%s%s?Hydrus-Client-API-Access-Key=%s&hash=%s', ...
            apiURL, searchEndpoint, apiKey, urlencode(fileHashes{potentialDuplicatePairs(i, 1)}));
        requestURL2 = sprintf('%s%s?Hydrus-Client-API-Access-Key=%s&hash=%s', ...
            apiURL, searchEndpoint, apiKey, urlencode(fileHashes{potentialDuplicatePairs(i, 2)}));
        try
            img1 = webread(requestURL1, options);
            img2 = webread(requestURL2, options);
            img1 = preprocessThumb(img1, P);
            img2 = preprocessThumb(img2, P);
            if P.useMSSSIM && exist('multissim', 'file')
                try
                    score = multissim(img1, img2);
                catch
                    score = ssim(img1, img2);
                end
            else
                score = ssim(img1, img2);
            end
            if score >= similarityThreshold
                duplicatesMatrix = [duplicatesMatrix; potentialDuplicatePairs(i, :)]; %#ok<AGROW>
            end
        catch ME
            if verbose
                fprintf('Error comparing thumbnails for pair %d: %s\n', i, ME.message);
            end
        end
    end
    if verbose
        fprintf('Completed comparison. Found %d duplicate pairs (threshold %.3f).\n', ...
            size(duplicatesMatrix, 1), similarityThreshold);
    end
end

% --- preprocessThumb ---
function img = preprocessThumb(img, P)
    if ~isempty(P.fixedResize)
        img = imresize(img, P.fixedResize);
    end
    if P.toGray && size(img,3) == 3
        img = rgb2gray(img);
    end
    if P.preBlurSigma > 0
        img = imgaussfilt(img, P.preBlurSigma);
    end
    if P.histEq
        img = histeq(img);
    end
end

% --- setFileRelationships ---
function report = setFileRelationships(api_key, api_url, file_hashes, pairs, opts)
    if nargin < 5 || isempty(opts), opts = struct(); end
    if isstring(file_hashes), file_hashes = cellstr(file_hashes); end
    if isempty(pairs), report = struct('attempted',0,'posted',0,'skipped',0,'failed',0,'fail_msgs',{{}}); return; end
    
    batch_size     = get_opt(opts,'batch_size',200);
    timeout_s      = get_opt(opts,'timeout_s',60);
    max_retries    = get_opt(opts,'max_retries',3);
    retry_pause_s  = get_opt(opts,'retry_pause_s',1.0);
    do_def_merge   = get_opt(opts,'do_default_merge',true);
    
    Pairs = unique(sort(pairs,2), 'rows', 'stable');
    endpoint = '/manage_file_relationships/set_file_relationships';
    url = [api_url, endpoint];
    
    baseOpts = weboptions('RequestMethod','POST','ContentType','json','Timeout',timeout_s, ...
                          'HeaderFields',{'Hydrus-Client-API-Access-Key', api_key});
    
    posted=0; skipped=0; failed=0; fail_msgs={};
    nPairs = size(Pairs,1);
    nBatches = ceil(nPairs / batch_size);
    
    for b = 1:nBatches
        sIdx = (b-1)*batch_size + 1;
        eIdx = min(b*batch_size, nPairs);
        chunk = Pairs(sIdx:eIdx, :);
        
        rels = repmat(struct('hash_a','', 'hash_b','', ...
                             'relationship', 0, ...
                             'do_default_content_merge', logical(do_def_merge)), size(chunk,1), 1);
        for k = 1:size(chunk,1)
            rels(k).hash_a = file_hashes{chunk(k,1)};
            rels(k).hash_b = file_hashes{chunk(k,2)};
        end
        
        body = struct('relationships', rels);
        [ok, msg] = post_json(url, body, baseOpts, max_retries, retry_pause_s);
        if ok
            posted = posted + numel(rels);
        else
            if contains(lower(msg),'400') || contains(lower(msg),'bad request')
                skipped = skipped + numel(rels);
            else
                failed = failed + numel(rels);
                fail_msgs{end+1} = msg; %#ok<AGROW>
            end
        end
    end
    
    report = struct('attempted', size(Pairs,1), 'posted', posted, 'skipped', skipped, 'failed', failed, 'fail_msgs', {fail_msgs});
end

function val = get_opt(s, f, d)
    if isfield(s,f) && ~isempty(s.(f)), val = s.(f); else, val = d; end
end

function [ok, errmsg] = post_json(url, body, options, max_retries, pause_s)
    ok = false; errmsg = '';
    for a = 1:max_retries
        try
            webwrite(url, body, options);
            ok = true; return
        catch err
            errmsg = err.message;
            if a < max_retries
                if contains(lower(errmsg),'400') || contains(lower(errmsg),'bad request')
                    break
                end
                pause(pause_s);
            end
        end
    end
end

% --- Tag fetch (modern + legacy Hydrus formats) ---
function tagLists = get_tags_for_hashes(hashes, api_key, api_URL, VERBOSE)
    endpoint   = '/get_files/file_metadata';
    batch_size = 256;   % Hydrus tends to use ~256 internally
    n          = numel(hashes);
    tagLists   = cell(1, n);

    opts = weboptions('RequestMethod','GET','Timeout',60);
    baseURL = sprintf([ ...
        '%s%s?' ...
        'only_return_basic_information=false&' ...
        'include_notes=false&' ...
        'include_services_object=true&' ...
        'Hydrus-Client-API-Access-Key=%s'], ...
        api_URL, endpoint, api_key);

    batches  = floor(n / batch_size);
    remainder = mod(n, batch_size);

    for b = 1:batches
        idx = (b-1)*batch_size + (1:batch_size);
        req = sprintf('%s&hashes=%s', baseURL, urlencode(jsonencode(hashes(idx))));
        R = webread(req, opts);
        tagLists(idx) = parseTagListsFromMetadata_FLEX(R);
        if VERBOSE
            fprintf('Fetched tags batch %d/%d (%d items)\n', b, batches + (remainder>0), batch_size);
        end
    end

    if remainder > 0
        idx = batches*batch_size + (1:remainder);
        req = sprintf('%s&hashes=%s', baseURL, urlencode(jsonencode(hashes(idx))));
        R = webread(req, opts);
        tagLists(idx) = parseTagListsFromMetadata_FLEX(R);
        if VERBOSE
            fprintf('Fetched tags remainder (%d items)\n', remainder);
        end
    end
end

function outCells = parseTagListsFromMetadata_FLEX(resp)
% Parse tags from /get_files/file_metadata for modern & legacy shapes.
% Robustly flattens nested structs/cells into one list of strings.
    m = resp.metadata;
    outCells = cell(1, numel(m));

    for k = 1:numel(m)
        allTags = {};                % collect as column for safe concatenation
        allTags = allTags(:);

        % Modern shape: per-service 'tags' object
        if isfield(m(k), 'tags') && isstruct(m(k).tags)
            svcKeys = fieldnames(m(k).tags);
            for s = 1:numel(svcKeys)
                svc = m(k).tags.(svcKeys{s});
                if isfield(svc, 'display_tags') && isstruct(svc.display_tags)
                    allTags = [allTags; flattenTagsAny(svc.display_tags)]; %#ok<AGROW>
                end
                if isfield(svc, 'storage_tags') && isstruct(svc.storage_tags)
                    allTags = [allTags; flattenTagsAny(svc.storage_tags)]; %#ok<AGROW>
                end
            end

        % Legacy: service_names_to_statuses_to_tags (name keyed)
        elseif isfield(m(k), 'service_names_to_statuses_to_tags') && isstruct(m(k).service_names_to_statuses_to_tags)
            map = m(k).service_names_to_statuses_to_tags;
            allTags = [allTags; flattenTagsAny(map)]; %#ok<AGROW>

        % Very old: flat 'tags' cell
        elseif isfield(m(k), 'tags') && iscell(m(k).tags)
            allTags = [allTags; flattenTagsAny(m(k).tags)]; %#ok<AGROW>
        end

        % Normalize to unique lowercase strings
        allTags = string(allTags(:));
        allTags = lower(allTags);
        allTags = allTags(~ismissing(allTags) & strlength(allTags) > 0);
        allTags = unique(allTags, 'stable');

        outCells{k} = cellstr(allTags);
    end
end

function out = flattenTagsAny(x)
% Recursively flatten structs / cells / strings into a column cellstr list.
    out = {};
    if isempty(x), return; end

    if isstruct(x)
        f = fieldnames(x);
        for i = 1:numel(f)
            out = [out; flattenTagsAny(x.(f{i}))]; %#ok<AGROW>
        end
        return
    end

    if iscell(x)
        for i = 1:numel(x)
            out = [out; flattenTagsAny(x{i})]; %#ok<AGROW>
        end
        return
    end

    if isstring(x)
        x = x(:);
        x = x(~ismissing(x));
        out = cellstr(x);
        return
    end

    if ischar(x)
        out = {x};
        return
    end
    % Other types ignored
end

function [tagsSorted, countsSorted] = countAllTags(tagLists)
% Count frequency of all tags across files (case-insensitive).
% Each tag counts once per file.
    tagCountMap = containers.Map('KeyType','char','ValueType','double');
    for k = 1:numel(tagLists)
        tags = tagLists{k};
        if isempty(tags), continue; end
        tags = unique(lower(string(tags(:))));
        for t = 1:numel(tags)
            tag = tags{t};
            if isKey(tagCountMap, tag)
                tagCountMap(tag) = tagCountMap(tag) + 1;
            else
                tagCountMap(tag) = 1;
            end
        end
    end
    tagsSorted = tagCountMap.keys;
    countsSorted = cell2mat(tagCountMap.values);
    [countsSorted, sortIdx] = sort(countsSorted, 'descend');
    tagsSorted = tagsSorted(sortIdx);
end

% --- aHash prefilter ---
function [pairsOut, stats] = prefilterPairsWithAHash(pairsIn, file_hash, api_key, api_URL, opts, VERBOSE)
% Prefilter candidate pairs using a tiny perceptual hash (aHash).
% Keeps pairs with Hamming distance <= opts.hammingMax (default 10).
    if nargin < 5 || isempty(opts), opts = struct(); end
    resize      = get_opt(opts,'resize',[64 64]);
    hammingMax  = get_opt(opts,'hammingMax',10);

    idx = unique(pairsIn(:));
    hashes = file_hash(idx);

    endpoint = '/get_files/thumbnail';
    wopt = weboptions('RequestMethod','GET','ContentType','image','Timeout',30);
    aHashes = containers.Map('KeyType','char','ValueType','uint64');

    for t = 1:numel(idx)
        h = hashes{t};
        url = sprintf('%s%s?Hydrus-Client-API-Access-Key=%s&hash=%s', api_URL, endpoint, api_key, urlencode(h));
        try
            im = webread(url, wopt);
            ah = computeAHash(im, resize);
            aHashes(h) = ah;
        catch ME
            if VERBOSE
                fprintf('aHash: failed thumb for idx %d: %s\n', idx(t), ME.message);
            end
        end
    end

    keep = false(size(pairsIn,1),1);
    droppedMissing = 0;
    for k = 1:size(pairsIn,1)
        i = pairsIn(k,1); j = pairsIn(k,2);
        hi = file_hash{i}; hj = file_hash{j};
        if ~isKey(aHashes,hi) || ~isKey(aHashes,hj)
            droppedMissing = droppedMissing + 1;
            continue
        end
        d = hamming64(aHashes(hi), aHashes(hj));
        keep(k) = (d <= hammingMax);
    end

    pairsOut = pairsIn(keep,:);
    stats = struct('missing', droppedMissing, 'kept', sum(keep), 'total', size(pairsIn,1));
end

function h = computeAHash(im, preResize)
% 8x8 average hash (aHash), returned as uint64 bitstring (MSB-first packing).
    if ~isempty(preResize)
        im = imresize(im, preResize);
    end
    if ndims(im) == 3 && size(im,3) == 3
        im = rgb2gray(im);
    end
    im = im2uint8(im);
    im = imresize(im, [8 8], 'nearest');

    m = mean(im(:), 'omitnan');
    bits = im > m;  % 8x8 logical

    h = uint64(0);
    for r = 1:8
        for c = 1:8
            if bits(r,c)
                shift = uint64(64 - ((r-1)*8 + c)); % 64..1
                h = bitor(h, bitshift(uint64(1), shift-1));
            end
        end
    end
end

function d = hamming64(a, b)
% Hamming distance between two uint64 bitstrings
    x = bitxor(a,b);
    d = 0;
    for k = 1:64
        d = d + bitand(x, uint64(1));
        x = bitshift(x, -1);
    end
    d = double(d);
end

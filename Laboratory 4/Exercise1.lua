local directory = arg[0]:match("^(.*)/") or "."
local genomePath = arg[1] or directory .. "/../Laboratory 3/GCF_002194405.1_ViralProj390412_genomic.fna"
local samplesPath = directory .. "/samples.txt"

local SAMPLES, MIN_LENGTH, MAX_LENGTH, MIN_OVERLAP = 20000, 100, 150, 10

local function readGenome(path)
    local lines = {}
    for line in io.lines(path) do
        if not line:match("^>") then table.insert(lines, line) end
    end
    return table.concat(lines):upper()
end

local function breakGenome(genome, path)
    local file = io.open(path, "w")
    local covered = {}
    for _ = 1, SAMPLES do
        local length = math.random(MIN_LENGTH, MAX_LENGTH)
        local start = math.random(1, #genome - length + 1)
        file:write(genome:sub(start, start + length - 1), "\n")
        for i = start, start + length - 1 do covered[i] = true end
    end
    file:close()
    local count = 0
    for i = 1, #genome do
        if covered[i] then count = count + 1 end
    end
    return count
end

local function readSamples(path)
    local samples = {}
    for line in io.lines(path) do table.insert(samples, line) end
    return samples
end

local function extend(contig, samples)
    local index = {}
    for i, sample in ipairs(samples) do
        local key = sample:sub(1, MIN_OVERLAP)
        index[key] = index[key] or {}
        table.insert(index[key], i)
    end

    local used = {}
    while true do
        local found, overlap
        for p = math.max(1, #contig - MAX_LENGTH + 2), #contig - MIN_OVERLAP + 1 do
            local tail = contig:sub(p)
            for _, i in ipairs(index[tail:sub(1, MIN_OVERLAP)] or {}) do
                local sample = samples[i]
                if not used[i] and #sample > #tail and sample:sub(1, #tail) == tail then
                    found, overlap = i, #tail
                    break
                end
            end
            if found then break end
        end
        if not found then return contig end
        used[found] = true
        contig = contig .. samples[found]:sub(overlap + 1)
    end
end

local function reverseAll(samples)
    local reversed = {}
    for i, sample in ipairs(samples) do reversed[i] = sample:reverse() end
    return reversed
end

math.randomseed(os.time())
local clock = os.clock()

local genome = readGenome(genomePath)
local covered = breakGenome(genome, samplesPath)
local samples = readSamples(samplesPath)

local right = extend(samples[1], samples)
local rebuilt = extend(right:reverse(), reverseAll(samples)):reverse()

print("Genome length:     " .. #genome)
print("Samples:           " .. #samples .. " (" .. MIN_LENGTH .. "-" .. MAX_LENGTH .. " letters, saved to samples.txt)")
print(string.format("Covered positions: %d / %d (%.2f%%)", covered, #genome, covered / #genome * 100))
print("Rebuilt length:    " .. #rebuilt)

if rebuilt == genome then
    print("Result:            rebuilt genome is identical to the original")
else
    local position = genome:find(rebuilt, 1, true)
    if position then
        print(string.format("Result:            matches the original from %d to %d (%d letters missing at the start, %d at the end)",
            position, position + #rebuilt - 1, position - 1, #genome - (position + #rebuilt - 1)))
    else
        print("Result:            rebuilt genome differs from the original")
    end
end
print(string.format("Time:              %.2f s", os.clock() - clock))

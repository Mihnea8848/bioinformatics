local directory = arg[0]:match("^(.*)/") or "."
local path = arg[1] or directory .. "/../Laboratory 3/GCF_002194405.1_ViralProj390412_genomic.fna"

local MIN_IR, MIN_TSD, MAX_TSD, MIN_CORE, MAX_LENGTH = 8, 4, 10, 5, 100
local COMPLEMENT = { A = "T", T = "A", C = "G", G = "C" }

local function reverseComplement(dna)
    return (dna:reverse():gsub(".", COMPLEMENT))
end

local function readFasta(path)
    local records, header, lines = {}, nil, {}
    local function save()
        if header then table.insert(records, { header = header, sequence = table.concat(lines):upper() }) end
    end
    for line in io.lines(path) do
        if line:match("^>") then
            save()
            header, lines = line:sub(2), {}
        else
            table.insert(lines, (line:gsub("%s", "")))
        end
    end
    save()
    return records
end

local function directRepeat(s, start, stop)
    for length = math.min(MAX_TSD, start - 1), MIN_TSD, -1 do
        local left = s:sub(start - length, start - 1)
        if left == s:sub(stop + 1, stop + length) then return left end
    end
end

local function findEnd(s, i)
    local target = reverseComplement(s:sub(i, i + MIN_IR - 1))
    for j = i + MIN_IR + MIN_CORE, math.min(#s - MIN_IR + 1, i + MAX_LENGTH) do
        local stop = j + MIN_IR - 1
        if s:sub(j, stop) == target then
            local tsd = directRepeat(s, i, stop)
            if tsd then return stop, tsd end
        end
    end
end

local function invertedRepeat(s, start, stop)
    local length = MIN_IR
    while COMPLEMENT[s:sub(start + length, start + length)] == s:sub(stop - length, stop - length) do length = length + 1 end
    return s:sub(start, start + length - 1)
end

local total = 0
for _, record in ipairs(readFasta(path)) do
    local s = record.sequence
    print(record.header)
    print("Length: " .. #s)
    print()
    print("#    Start     End  Inverted repeat           Direct repeat")

    local i, found = 1, 0
    while i <= #s - MIN_IR do
        local stop, tsd = findEnd(s, i)
        if stop then
            found = found + 1
            local ir = invertedRepeat(s, i, stop)
            print(string.format("%-3d %6d  %6d  %s/%s  %s", found, i, stop, ir, reverseComplement(ir), tsd))
            i = stop + 1
        else
            i = i + 1
        end
    end
    print()
    total = total + found
end

print("Transposable elements found: " .. total)

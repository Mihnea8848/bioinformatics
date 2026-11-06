local directory = arg[0]:match("^(.*)/") or "."
local path = arg[1] or directory .. "/sequence.txt"

local MIN_IR, MIN_TSD, MAX_TSD, MIN_CORE, MAX_LENGTH = 8, 4, 10, 5, 100
local COMPLEMENT = { A = "T", T = "A", C = "G", G = "C" }

local function reverseComplement(dna)
    return (dna:reverse():gsub(".", COMPLEMENT))
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

local s = io.open(path):read("a"):upper():gsub("%s", "")
print("Sequence length: " .. #s)
print()
print("#  Start  End  Inverted repeat  Direct repeat")

local i, found = 1, 0
while i <= #s - MIN_IR do
    local stop, tsd = findEnd(s, i)
    if stop then
        found = found + 1
        local ir = invertedRepeat(s, i, stop)
        print(string.format("%d  %5d  %3d  %-8s/%-8s %s", found, i, stop, ir, reverseComplement(ir), tsd))
        i = stop + 1
    else
        i = i + 1
    end
end

print()
print("Transposable elements found: " .. found)

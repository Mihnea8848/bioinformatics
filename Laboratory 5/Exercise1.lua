local directory = arg[0]:match("^(.*)/") or "."
local path = directory .. "/sequence.txt"

local IR_LENGTH, TSD_LENGTH, MIN_CORE, MAX_CORE = 8, 5, 10, 30
local BASES = "ACGT"
local COMPLEMENT = { A = "T", T = "A", C = "G", G = "C" }

local function randomDNA(length)
    local bases = {}
    for i = 1, length do
        local k = math.random(#BASES)
        bases[i] = BASES:sub(k, k)
    end
    return table.concat(bases)
end

local function reverseComplement(dna)
    return (dna:reverse():gsub(".", COMPLEMENT))
end

math.randomseed(os.time())

local transposons, inserted = {}, 0
for i = 1, math.random(3, 4) do
    local ir = randomDNA(IR_LENGTH)
    transposons[i] = ir .. randomDNA(math.random(MIN_CORE, MAX_CORE)) .. reverseComplement(ir)
    inserted = inserted + #transposons[i] + TSD_LENGTH
end

local host = randomDNA(math.random(math.max(200, inserted + 100), 400) - inserted)
local segment = #host // #transposons

local sequence, previous, elements = "", 0, {}
for i, transposon in ipairs(transposons) do
    local site = math.random((i - 1) * segment + TSD_LENGTH + 2, i * segment - 2)
    local tsd = host:sub(site - TSD_LENGTH + 1, site)
    sequence = sequence .. host:sub(previous + 1, site)
    elements[i] = { start = #sequence + 1, stop = #sequence + #transposon, ir = transposon:sub(1, IR_LENGTH), tsd = tsd }
    sequence = sequence .. transposon .. tsd
    previous = site
end
sequence = sequence .. host:sub(previous + 1)

local file = io.open(path, "w")
file:write(sequence, "\n")
file:close()

print("Sequence (" .. #sequence .. " bases, saved to sequence.txt):")
print(sequence)
print()
print("#  Start  End  Inverted repeat  Direct repeat")
for i, element in ipairs(elements) do
    print(string.format("%d  %5d  %3d  %-8s/%-8s %s", i, element.start, element.stop, element.ir, reverseComplement(element.ir), element.tsd))
end

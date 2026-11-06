local directory = arg[0]:match("^(.*)/") or "."
local path = arg[1] or directory .. "/GCF_002194405.1_ViralProj390412_genomic.fna"

local bases = "UCAG"
local counts = {}
for first in bases:gmatch(".") do
    for second in bases:gmatch(".") do
        for third in bases:gmatch(".") do
            counts[first .. second .. third] = 0
        end
    end
end

local lines = {}
for line in io.lines(path) do
    if not line:match("^>") then table.insert(lines, line) end
end
local sequence = table.concat(lines):upper():gsub("%s", ""):gsub("T", "U")

local total = 0
for i = 1, #sequence - 2 do
    local codon = sequence:sub(i, i + 2)
    if counts[codon] then
        counts[codon] = counts[codon] + 1
        total = total + 1
    end
end

local codons = {}
for codon in pairs(counts) do table.insert(codons, codon) end
table.sort(codons, function(a, b) return counts[a] > counts[b] or (counts[a] == counts[b] and a < b) end)

print("File:   " .. path)
print("Length: " .. #sequence)
print("Codons: " .. total)
print()
print("Rank  Codon  Count  Frequency")
for rank = 1, 10 do
    local codon = codons[rank]
    print(string.format("%4d  %-5s  %5d  %8.2f%%", rank, codon, counts[codon], total > 0 and counts[codon] / total * 100 or 0))
end

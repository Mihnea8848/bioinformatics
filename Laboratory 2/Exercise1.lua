local S = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

local bases = "UCAG"
local aminoAcids = {
    "Phe", "Phe", "Leu", "Leu", "Ser", "Ser", "Ser", "Ser", "Tyr", "Tyr", "Stop", "Stop", "Cys", "Cys", "Stop", "Trp",
    "Leu", "Leu", "Leu", "Leu", "Pro", "Pro", "Pro", "Pro", "His", "His", "Gln", "Gln", "Arg", "Arg", "Arg", "Arg",
    "Ile", "Ile", "Ile", "Met", "Thr", "Thr", "Thr", "Thr", "Asn", "Asn", "Lys", "Lys", "Ser", "Ser", "Arg", "Arg",
    "Val", "Val", "Val", "Val", "Ala", "Ala", "Ala", "Ala", "Asp", "Asp", "Glu", "Glu", "Gly", "Gly", "Gly", "Gly",
}

local geneticCode = {}
local i = 0
for first in bases:gmatch(".") do
    for second in bases:gmatch(".") do
        for third in bases:gmatch(".") do
            i = i + 1
            geneticCode[first .. second .. third] = aminoAcids[i]
        end
    end
end

local function translate(dna)
    local rna = dna:upper():gsub("T", "U")
    local protein = {}
    for codon in rna:gmatch("...") do
        local aminoAcid = geneticCode[codon]
        if aminoAcid == "Stop" then break end
        table.insert(protein, aminoAcid)
    end
    return table.concat(protein, "-")
end

print("DNA:     " .. S)
print("Protein: " .. translate(S))

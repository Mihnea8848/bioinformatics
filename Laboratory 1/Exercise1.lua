-- Detect the frequency of each letter of the DNA alphabet in a sequence S.

local S = "ATCGACGCTAACGTA"

local function countLetters(sequence)
    local counts = {}
    for letter in sequence:gmatch(".") do
        counts[letter] = (counts[letter] or 0) + 1
    end
    return counts
end

local function getAlphabet(counts)
    local alphabet = {}
    for letter in pairs(counts) do
        table.insert(alphabet, letter)
    end
    table.sort(alphabet)
    return alphabet
end

local counts = countLetters(S)
local alphabet = getAlphabet(counts)
local total = #S

print("Sequence: " .. S)
print("Length:   " .. total)
print("Alphabet: " .. table.concat(alphabet, ", "))
print()
print("Letter  Count  Frequency")
for _, letter in ipairs(alphabet) do
    local count = counts[letter]
    print(string.format("%-6s  %5d  %8.2f%%", letter, count, count / total * 100))
end

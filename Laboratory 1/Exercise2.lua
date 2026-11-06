-- Measure the frequency of combinations of two and three letters in a sequence S, display on a chart

local S = "ATCGACGCTAACGTA"

local function countCombinations(sequence, k)
    local counts = {}
    local total = #sequence - k + 1
    for i = 1, total do
        local combination = sequence:sub(i, i + k - 1)
        counts[combination] = (counts[combination] or 0) + 1
    end
    return counts, total
end

local function showFrequencies(title, counts, total)
    local keys = {}
    for key in pairs(counts) do
        table.insert(keys, key)
    end
    table.sort(keys)

    print(title)
    for _, key in ipairs(keys) do
        print(string.format("%-4s %6.2f%%", key, counts[key] / total * 100))
    end
    print()

    -- gnuplot?
    local gnuplot = io.popen("gnuplot -persist", "w")
    gnuplot:write("set title '" .. title .. "'\n")
    gnuplot:write("set ylabel 'Frequency (%)'\n")
    gnuplot:write("set style fill solid\n")
    gnuplot:write("set boxwidth 0.7\n")
    gnuplot:write("set yrange [0:*]\n")
    gnuplot:write("plot '-' using 2:xtic(1) with boxes notitle\n")
    for _, key in ipairs(keys) do
        gnuplot:write(key .. " " .. counts[key] / total * 100 .. "\n")
    end
    gnuplot:write("e\n")
    gnuplot:close()
end

showFrequencies("Letters", countCombinations(S, 1))
showFrequencies("Combinations of two letters", countCombinations(S, 2))
showFrequencies("Combinations of three letters", countCombinations(S, 3))

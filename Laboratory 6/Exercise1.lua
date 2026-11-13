local directory = arg[0]:match("^(.*)/") or "."
local path = arg[1] or directory .. "/../Laboratory 3/GCF_002194405.1_ViralProj390412_genomic.fna"

local ENZYMES = {
    EcoRI = { site = "GAATTC", cut = 1 },
    BamHI = { site = "GGATCC", cut = 1 },
    HindIII = { site = "AAGCTT", cut = 1 },
    TaqI = { site = "TCGA", cut = 1 },
    HaeIII = { site = "GGCC", cut = 2 },
}
local SELECTED = #arg > 1 and { table.unpack(arg, 2) } or { "EcoRI", "BamHI", "HindIII" }
local LADDER = { 100, 200, 500, 1000, 2000, 3000, 5000, 10000, 20000 }
local ROWS, LANE = 30, 10

local function readGenome(path)
    local lines = {}
    for line in io.lines(path) do
        if not line:match("^>") then table.insert(lines, line) end
    end
    return table.concat(lines):upper()
end

local function cleave(dna, enzyme)
    local cuts, start = {}, dna:find(enzyme.site, 1, true)
    while start do
        table.insert(cuts, start + enzyme.cut - 1)
        start = dna:find(enzyme.site, start + 1, true)
    end
    return cuts
end

local function fragments(cuts, length)
    local result, previous = {}, 0
    for _, cut in ipairs(cuts) do
        table.insert(result, cut - previous)
        previous = cut
    end
    table.insert(result, length - previous)
    return result
end

local function gel(lanes, length)
    local top, bottom = math.log(math.max(length, LADDER[#LADDER])), math.log(LADDER[1])
    local function row(size)
        local r = math.floor((top - math.log(size)) / (top - bottom) * (ROWS - 1) + 0.5) + 1
        return math.max(1, math.min(ROWS, r))
    end

    local labels, bands = {}, {}
    for _, size in ipairs(LADDER) do labels[row(size)] = size end
    for i, lane in ipairs(lanes) do
        bands[i] = {}
        for _, size in ipairs(lane.fragments) do
            if size >= LADDER[1] then bands[i][row(size)] = (bands[i][row(size)] or 0) + 1 end
        end
    end

    local header, wells = string.rep(" ", 8), string.rep(" ", 8)
    for _, lane in ipairs(lanes) do
        header = header .. string.format("%-" .. LANE .. "s", lane.name)
        wells = wells .. "[" .. string.rep(" ", LANE - 4) .. "]  "
    end
    print(header)
    print(wells)
    for r = 1, ROWS do
        local line = string.format("%7s ", labels[r] or "")
        for i in ipairs(lanes) do
            local count = bands[i][r]
            local band = not count and " " or count == 1 and "-" or count == 2 and "=" or "#"
            line = line .. string.rep(band, LANE - 2) .. "  "
        end
        print(line)
    end
end

local genome = readGenome(path)
print("Genome length: " .. #genome)
print()

local lanes, all = { { name = "Ladder", fragments = LADDER } }, {}
for _, name in ipairs(SELECTED) do
    local enzyme = assert(ENZYMES[name], "unknown enzyme " .. name)
    local cuts = cleave(genome, enzyme)
    local pieces = fragments(cuts, #genome)
    print(string.format("%s (%s^%s): %d cleavages", name, enzyme.site:sub(1, enzyme.cut), enzyme.site:sub(enzyme.cut + 1), #cuts))
    print("  Positions: " .. table.concat(cuts, ", "))
    print("  Fragments: " .. table.concat(pieces, ", "))
    print()
    table.insert(lanes, { name = name, fragments = pieces })
    for _, cut in ipairs(cuts) do all[cut] = true end
end

local combined = {}
for cut in pairs(all) do table.insert(combined, cut) end
table.sort(combined)
table.insert(lanes, { name = "All", fragments = fragments(combined, #genome) })
print(string.format("All enzymes: %d cleavages, %d fragments", #combined, #combined + 1))
print()

gel(lanes, #genome)

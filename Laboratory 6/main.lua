local ENZYMES = {
    { name = "EcoRI", site = "GAATTC", cut = 1, color = { 0.84, 0.37, 0.00 }, enabled = true },
    { name = "BamHI", site = "GGATCC", cut = 1, color = { 0.00, 0.45, 0.70 }, enabled = true },
    { name = "HindIII", site = "AAGCTT", cut = 1, color = { 0.00, 0.62, 0.45 }, enabled = true },
    { name = "TaqI", site = "TCGA", cut = 1, color = { 0.90, 0.62, 0.00 }, enabled = false },
    { name = "HaeIII", site = "GGCC", cut = 2, color = { 0.80, 0.47, 0.65 }, enabled = false },
}
local COMPLEMENT = { A = "T", T = "A", C = "G", G = "C" }
local LADDER = { 100, 200, 500, 1000, 2000, 3000, 5000, 10000, 20000 }
local DEFAULT_GENOME = "/../Laboratory 3/GCF_002194405.1_ViralProj390412_genomic.fna"
local FONT_DIR = "/usr/share/fonts/liberation/"
local ROW, SMALLEST, LANE_WIDTH = 20, 80, 62
local INK, GREY, RULE, SHADE = { 0.1, 0.1, 0.1 }, { 0.42, 0.42, 0.42 }, { 0.8, 0.8, 0.8 }, { 0.94, 0.94, 0.94 }
local WELL, ERROR = { 0.04, 0.04, 0.05 }, { 0.75, 0.15, 0.1 }
local SELECTION, HIGHLIGHT = { 1, 0.78, 0.2 }, { 1, 0.93, 0.62 }

local fonts, textures, genome, lanes, targets, combined = {}, {}, nil, {}, {}, nil
local focus, scroll, progress, message, hover, list = 2, 0, 0, nil, nil, nil
local selected, hoveredBand = nil, nil

local function commas(n)
    return (tostring(math.floor(n)):reverse():gsub("(%d%d%d)", "%1,"):reverse():gsub("^,", ""))
end

local function plural(n, word)
    return commas(n) .. " " .. word .. (n == 1 and "" or "s")
end

local function text(str, x, y, font, color, align, width)
    love.graphics.setFont(font)
    love.graphics.setColor(color)
    love.graphics.printf(str, x, y, width or 10000, align or "left")
end

local function fit(str, font, width)
    if font:getWidth(str) <= width then return str end
    while #str > 0 and font:getWidth(str .. "...") > width do str = str:sub(1, -2) end
    return str .. "..."
end

local function rule(x, y, w, weight, color)
    love.graphics.setColor(color or INK)
    love.graphics.setLineWidth(weight)
    love.graphics.line(x, y, x + w, y)
end

local function inside(x, y, w, h)
    local mx, my = love.mouse.getPosition()
    return mx >= x and mx <= x + w and my >= y and my <= y + h
end

local function loadFont(name, size)
    local file = io.open(FONT_DIR .. name, "rb")
    if not file then return love.graphics.newFont(size) end
    local data = love.filesystem.newFileData(file:read("*a"), name)
    file:close()
    return love.graphics.newFont(love.font.newTrueTypeRasterizer(data, size))
end

local function readFasta(path)
    local file = io.open(path)
    if not file then return nil, "Cannot open " .. path:match("[^/]+$") end
    local title, parts = nil, {}
    for line in file:lines() do
        if line:sub(1, 1) == ">" then
            title = title or line:sub(2)
        else
            parts[#parts + 1] = line:gsub("%s", "")
        end
    end
    file:close()
    local dna = table.concat(parts):upper()
    if #dna == 0 then return nil, "No DNA sequence in " .. path:match("[^/]+$") end
    return title or path:match("[^/]+$"), dna
end

local function cleave(dna, enzyme)
    local cuts, start = {}, dna:find(enzyme.site, 1, true)
    while start do
        cuts[#cuts + 1] = start + enzyme.cut - 1
        start = dna:find(enzyme.site, start + 1, true)
    end
    return cuts
end

local function overhang(enzyme)
    local cut, other = enzyme.cut, #enzyme.site - enzyme.cut
    if cut == other then return "blunt" end
    if cut < other then return "5' " .. enzyme.site:sub(cut + 1, other) end
    return "3' " .. enzyme.site:sub(other + 1, cut)
end

local function group(pieces)
    local bySize, bands = {}, {}
    for index, piece in ipairs(pieces) do
        local band = bySize[piece.size]
        if not band then
            band = { size = piece.size, count = 0, first = index }
            bySize[piece.size], bands[#bands + 1] = band, band
        end
        band.count = band.count + 1
    end
    table.sort(bands, function(a, b) return a.size > b.size end)
    return bands
end

local function makeTextures()
    local band = love.image.newImageData(64, 32)
    band:mapPixel(function(x, y)
        local dy, dx = (y + 0.5 - 16) / 16, math.abs((x + 0.5 - 32) / 32)
        local edge = dx < 0.82 and 1 or math.max(0, 1 - (dx - 0.82) / 0.18)
        return 1, 1, 1, math.exp(-dy * dy * 7) * edge
    end)
    textures.band = love.graphics.newImage(band)

    local noise = love.image.newImageData(256, 256)
    noise:mapPixel(function()
        local value = 0.13 + love.math.random() * 0.035
        return value, value, value + 0.005, 1
    end)
    textures.noise = love.graphics.newImage(noise)
    textures.noise:setWrap("repeat", "repeat")
    textures.quad = love.graphics.newQuad(0, 0, 1, 1, 256, 256)
end

local function round3(value)
    local power = 10 ^ (math.floor(math.log10(value)) - 2)
    return math.floor(value / power + 0.5) * power
end

local function makeLane(name, short, cuts, ownerOf)
    local pieces, previous = {}, 0
    for _, cut in ipairs(cuts) do
        pieces[#pieces + 1] = { start = previous + 1, stop = cut, size = cut - previous,
            left = previous > 0 and ownerOf(previous) or "terminus", right = ownerOf(cut) }
        previous = cut
    end
    pieces[#pieces + 1] = { start = previous + 1, stop = #genome.dna, size = #genome.dna - previous,
        left = previous > 0 and ownerOf(previous) or "terminus", right = "terminus" }
    return { name = name, short = short, pieces = pieces, bands = group(pieces) }
end

local function rebuild()
    local ladder = {}
    for i, size in ipairs(LADDER) do ladder[i] = { size = size } end
    lanes = { { name = "DNA size marker", short = "size marker", pieces = ladder, bands = group(ladder), ladder = true } }

    local owners, cuts, names = {}, {}, {}
    for _, enzyme in ipairs(ENZYMES) do
        if enzyme.enabled then
            lanes[#lanes + 1] = makeLane(enzyme.name, enzyme.name, enzyme.cuts, function() return enzyme.name end)
            names[#names + 1] = enzyme.name
            for _, cut in ipairs(enzyme.cuts) do owners[cut] = owners[cut] or enzyme end
        end
    end
    combined = nil
    if #names > 1 then
        for cut in pairs(owners) do cuts[#cuts + 1] = cut end
        table.sort(cuts)
        combined = { name = "All", cuts = cuts, owners = owners }
        lanes[#lanes + 1] = makeLane(table.concat(names, " + "), "combined digest", cuts, function(cut) return owners[cut].name end)
    end
    for i, lane in ipairs(lanes) do lane.tag = i == 1 and "M" or tostring(i - 1) end
    focus, scroll, progress, selected = math.min(focus, #lanes), 0, 0, nil
end

local function loadGenome(path)
    local title, dna = readFasta(path)
    if not title then
        message = dna
        return
    end
    local gc = select(2, dna:gsub("[GC]", "")) / #dna * 100
    genome, message = { title = title, dna = dna, gc = gc }, nil
    for _, enzyme in ipairs(ENZYMES) do enzyme.cuts = cleave(dna, enzyme) end
    focus = 2
    rebuild()
end

local function heading(letter, title, x, y, w, note)
    text(letter, x, y, fonts.heading, INK)
    text(title, x + 22, y, fonts.heading, INK)
    if note then text(note, x, y + 3, fonts.small, GREY, "right", w) end
end

local function checkbox(x, y, checked)
    love.graphics.setColor(INK)
    love.graphics.setLineWidth(1)
    love.graphics.rectangle("line", x + 0.5, y + 0.5, 12, 12)
    if checked then
        love.graphics.setLineWidth(2)
        love.graphics.line(x + 3, y + 6.5, x + 5.5, y + 9.5, x + 10, y + 3)
    end
end

local function drawHeader(w)
    text("Restriction Digest Analysis", 24, 16, fonts.title, INK)
    local info = genome and string.format("Sequence: %s     Length: %s bp (linear)     GC content: %.1f%%",
        genome.title, commas(#genome.dna), genome.gc) or "No sequence loaded"
    text(fit(info, fonts.body, w * 0.66), 24, 50, fonts.body, GREY)
    text(message or "Drag a FASTA file onto the window to analyse another sequence", w * 0.4, 52, fonts.small,
        message and ERROR or GREY, "right", w * 0.6 - 24)
    rule(24, 76, w - 48, 1.5)
end

local function drawEnzymes(x, y, w)
    heading("A", "Restriction enzymes", x, y, w, "select the enzymes to digest with")
    local top, rowHeight, cw = y + 30, 44, fonts.mono:getWidth("A")
    rule(x, top, w, 1.2)
    for _, column in ipairs({ { "Enzyme", 22 }, { "Recognition site", 136 }, { "Ends", 286 } }) do
        text(column[1], x + column[2], top + 5, fonts.smallBold, INK)
    end
    text("Sites", x, top + 5, fonts.smallBold, INK, "right", w - 4)
    rule(x, top + 24, w, 0.8)

    for i, enzyme in ipairs(ENZYMES) do
        local ry = top + 26 + (i - 1) * rowHeight
        if inside(x, ry, w, rowHeight) then
            love.graphics.setColor(SHADE)
            love.graphics.rectangle("fill", x, ry, w, rowHeight)
        end
        checkbox(x + 2, ry + 15, enzyme.enabled)
        love.graphics.setColor(enzyme.color)
        love.graphics.rectangle("fill", x + 22, ry + 17, 8, 8)
        text(enzyme.name, x + 36, ry + 13, fonts.body, INK)

        local seqX, length = x + 156, #enzyme.site
        local bottom = enzyme.site:gsub(".", COMPLEMENT)
        text("5'", seqX - 18, ry + 6, fonts.small, GREY)
        text(enzyme.site, seqX, ry + 4, fonts.mono, INK)
        text("3'", seqX + length * cw + 6, ry + 6, fonts.small, GREY)
        text("3'", seqX - 18, ry + 24, fonts.small, GREY)
        text(bottom, seqX, ry + 22, fonts.mono, INK)
        text("5'", seqX + length * cw + 6, ry + 24, fonts.small, GREY)
        local top, bottomCut = seqX + enzyme.cut * cw, seqX + (length - enzyme.cut) * cw
        love.graphics.setColor(enzyme.color)
        love.graphics.setLineWidth(1.5)
        love.graphics.line(top, ry + 4, top, ry + 21, bottomCut, ry + 21, bottomCut, ry + 39)

        text(overhang(enzyme), x + 286, ry + 13, fonts.mono, INK)
        text(commas(#enzyme.cuts), x, ry + 13, fonts.mono, INK, "right", w - 4)

        targets[#targets + 1] = { x = x, y = ry, w = w, h = rowHeight, action = function()
            enzyme.enabled = not enzyme.enabled
            rebuild()
        end }
    end
    local bottom = top + 26 + #ENZYMES * rowHeight + 2
    rule(x, bottom, w, 1.2)
    return bottom
end

local function drawFragments(x, y, w, h)
    local lane = lanes[focus]
    heading("B", "Fragments in lane " .. lane.tag, x, y, w, plural(#lane.pieces, "fragment"))
    text(fit(lane.name, fonts.small, w - 22), x + 22, y + 20, fonts.small, GREY)

    local top = y + 40
    local columns = {
        { "No.", 0, 30, "right" }, { "Start", 30, 96, "right" }, { "End", 96, 162, "right" },
        { "Length (bp)", 162, 246, "right" }, { "Left end", 262, 344, "left" }, { "Right end", 344, w, "left" },
    }
    local function cell(value, column, cy, font, color)
        local width = column[3] - column[2] - (column[4] == "right" and 8 or 0)
        text(value, x + column[2], cy, font, color, column[4], width)
    end
    rule(x, top, w, 1.2)
    for _, column in ipairs(columns) do cell(column[1], column, top + 5, fonts.smallBold, INK) end
    rule(x, top + 24, w, 0.8)

    local listTop, height = top + 26, h - 26 - 42
    list = { x = x, y = listTop, w = w, h = height, max = math.max(0, #lane.pieces * ROW - height) }
    scroll = math.min(scroll, list.max)
    love.graphics.setScissor(x, listTop, w, height)
    local first = math.floor(scroll / ROW) + 1
    for i = first, math.min(#lane.pieces, first + math.ceil(height / ROW)) do
        local piece, ry = lane.pieces[i], listTop + (i - 1) * ROW - scroll
        local chosen = selected and selected.lane == focus and selected.size == piece.size
        if chosen or i % 2 == 0 or inside(x, ry, w, ROW) then
            love.graphics.setColor(chosen and HIGHLIGHT or SHADE)
            love.graphics.rectangle("fill", x, ry, w, ROW)
        end
        if ry + ROW > listTop and ry < listTop + height then
            targets[#targets + 1] = { x = x, y = math.max(ry, listTop), w = w, h = ROW, action = function()
                selected = { lane = focus, size = piece.size }
            end }
        end
        local values = {
            i, piece.start and commas(piece.start) or "-", piece.stop and commas(piece.stop) or "-",
            commas(piece.size), piece.left or "-", piece.right or "-",
        }
        for c, column in ipairs(columns) do
            cell(values[c], column, ry + 3, fonts.monoSmall, (c == 1 or values[c] == "terminus") and GREY or INK)
        end
    end
    love.graphics.setScissor()
    rule(x, listTop + height + 2, w, 1.2)
end

local function drawLegend(x, y, w)
    rule(x, y, w, 1.2)
    text("Lane", x, y + 5, fonts.smallBold, INK)
    text("Digest", x + 40, y + 5, fonts.smallBold, INK)
    text("Bands", x, y + 5, fonts.smallBold, INK, "right", w - 4)
    rule(x, y + 24, w, 0.8)
    for i, lane in ipairs(lanes) do
        local ry = y + 26 + (i - 1) * 22
        if i == focus or inside(x, ry, w, 22) then
            love.graphics.setColor(SHADE)
            love.graphics.rectangle("fill", x, ry, w, 22)
        end
        local font = i == focus and fonts.bodyBold or fonts.body
        text(lane.tag, x + 4, ry + 4, font, INK)
        text(fit(lane.name, font, w - 90), x + 40, ry + 4, font, INK)
        text(commas(#lane.bands), x, ry + 4, fonts.monoSmall, INK, "right", w - 4)
        targets[#targets + 1] = { x = x, y = ry, w = w, h = 22, action = function() focus, scroll = i, 0 end }
    end
    local bottom = y + 26 + #lanes * 22 + 2
    rule(x, bottom, w, 1.2)
    text("Migration distance is modelled as linear in log10(fragment length). Band intensity is proportional to DNA mass. "
        .. "Fragments shorter than " .. SMALLEST .. " bp have run off the gel. Sizes of the selected lane are labelled on the right of the gel.",
        x, bottom + 10, fonts.caption, GREY, "left", w)
end

local function dashed(x1, x2, y)
    for dx = x1, x2, 8 do love.graphics.line(dx, y, math.min(dx + 4, x2), y) end
end

local function drawAnnotations(lane, x, top, bottom)
    local labels, spacing = {}, 14
    for _, band in ipairs(lane.bands) do
        if band.y then labels[#labels + 1] = { y = band.y, ly = band.y, band = band } end
    end
    if #labels == 0 then return end
    if #labels * spacing > bottom - top + 24 then
        text(plural(#labels, "band") .. "\ntoo many\nto label", x + 8, top, fonts.small, GREY)
        return
    end
    for k = 2, #labels do labels[k].ly = math.max(labels[k].y, labels[k - 1].ly + spacing) end
    if labels[#labels].ly > bottom + 8 then
        labels[#labels].ly = bottom + 8
        for k = #labels - 1, 1, -1 do labels[k].ly = math.min(labels[k].ly, labels[k + 1].ly - spacing) end
    end
    love.graphics.setLineWidth(1)
    for _, label in ipairs(labels) do
        local chosen = selected and selected.lane == focus and selected.size == label.band.size
        love.graphics.setColor(INK)
        love.graphics.polygon("fill", x + 1, label.y, x + 6, label.y - 3, x + 6, label.y + 3)
        love.graphics.line(x + 6, label.y, x + 14, label.ly)
        local value = commas(label.band.size) .. (label.band.count > 1 and " x" .. label.band.count or "")
        if chosen then
            love.graphics.setColor(HIGHLIGHT)
            love.graphics.rectangle("fill", x + 15, label.ly - 7, fonts.monoSmall:getWidth(value) + 6, 14)
        end
        text(value, x + 18, label.ly - 7, fonts.monoSmall, INK)
    end
end

local function drawGel(x, y, w, h)
    heading("C", "Simulated agarose gel electrophoresis", x, y, w, "hover to measure, click a band to select it")
    local laneWidth = math.min(LANE_WIDTH, (w * 0.55 - 80) / #lanes)
    local gelX, gelTop, gelH = x + 62, y + 58, h - 66
    local gelW = #lanes * laneWidth + 12
    local wellY, top, bottom = gelTop + 12, gelTop + 30, gelTop + gelH - 14
    local bandHeight = math.max(10, math.min(16, gelH / 26))

    local logMax, logMin = math.log(math.max(#genome.dna, LADDER[#LADDER]) * 1.15), math.log(SMALLEST)
    local function position(size)
        return top + (logMax - math.log(size)) / (logMax - logMin) * (bottom - top)
    end
    local function sizeAt(py)
        return math.exp(logMax - (py - top) / (bottom - top) * (logMax - logMin))
    end
    local maxMass = 1
    for i = 2, #lanes do
        for _, band in ipairs(lanes[i].bands) do maxMass = math.max(maxMass, band.size * band.count) end
    end

    love.graphics.setColor(1, 1, 1)
    textures.quad:setViewport(0, 0, gelW, gelH, 256, 256)
    love.graphics.draw(textures.noise, textures.quad, gelX, gelTop)
    text("bp", x, gelTop - 18, fonts.small, GREY, "right", 50)
    love.graphics.setLineWidth(1)
    for _, size in ipairs(LADDER) do
        local ty = position(size)
        text(commas(size), x, ty - 7, fonts.monoSmall, INK, "right", 50)
        love.graphics.setColor(INK)
        love.graphics.line(gelX - 8, ty, gelX - 2, ty)
    end

    local eased, mx, my = 1 - (1 - progress) ^ 3, love.mouse.getPosition()
    hoveredBand = nil
    for i, lane in ipairs(lanes) do
        local lx = gelX + 6 + (i - 1) * laneWidth
        local cx, bw = lx + laneWidth / 2, laneWidth * 0.74
        local font = i == focus and fonts.bodyBold or fonts.body
        text(lane.tag, lx, gelTop - 20, font, INK, "center", laneWidth)
        if i == focus then rule(cx - 7, gelTop - 4, 14, 1.5) end
        love.graphics.setColor(WELL)
        love.graphics.rectangle("fill", cx - bw / 2, wellY - 4, bw, 8)

        love.graphics.setBlendMode("add")
        if not lane.ladder then
            love.graphics.setColor(1, 1, 1, 0.02 * eased)
            love.graphics.rectangle("fill", cx - bw / 2 + 2, wellY + 4, bw - 4, (bottom - wellY) * eased)
        end
        for _, band in ipairs(lane.bands) do
            band.y = nil
            if band.size >= SMALLEST then
                band.y = wellY + (position(band.size) - wellY) * eased
                local intensity = lane.ladder and 0.75
                    or 0.35 + 0.65 * math.max(0, math.log(band.size * band.count / SMALLEST) / math.log(maxMass / SMALLEST))
                love.graphics.setColor(1, 1, 1, intensity)
                love.graphics.draw(textures.band, cx - bw / 2, band.y - bandHeight / 2, 0, bw / 64, bandHeight / 32)
                if progress == 1 and math.abs(mx - cx) <= bw / 2 and math.abs(my - band.y) <= bandHeight / 3 then
                    hoveredBand = { lane = i, band = band }
                end
            end
        end
        love.graphics.setBlendMode("alpha")

        for _, band in ipairs(lane.bands) do
            local chosen = selected and selected.lane == i and selected.size == band.size
            local over = hoveredBand and hoveredBand.band == band
            if band.y and (chosen or over) then
                love.graphics.setColor(chosen and SELECTION or { 1, 1, 1, 0.7 })
                love.graphics.setLineWidth(chosen and 2 or 1)
                love.graphics.rectangle("line", cx - bw / 2 - 3, band.y - bandHeight / 2, bw + 6, bandHeight)
            end
        end

        targets[#targets + 1] = { x = lx, y = gelTop - 22, w = laneWidth, h = gelH + 22, action = function()
            if hoveredBand and hoveredBand.lane == i then
                focus, selected = i, { lane = i, size = hoveredBand.band.size }
                scroll = (hoveredBand.band.first - 1) * ROW
            else
                focus, scroll, selected = i, 0, nil
            end
        end }
    end

    if hoveredBand then
        local lane, band = lanes[hoveredBand.lane], hoveredBand.band
        hover = { string.format("Lane %s (%s)  %s bp%s", lane.tag, lane.short, commas(band.size),
            band.count > 1 and ", " .. band.count .. " fragments" or ""), mx, my }
    elseif progress == 1 and inside(gelX, top - 10, gelW, bottom - top + 20) then
        love.graphics.setColor(1, 1, 1, 0.5)
        love.graphics.setLineWidth(1)
        dashed(gelX, gelX + gelW, my)
        local label = "~" .. commas(round3(sizeAt(my)))
        local lw = fonts.monoSmall:getWidth(label) + 8
        love.graphics.setColor(1, 1, 1)
        love.graphics.rectangle("fill", gelX - lw - 4, my - 9, lw, 18)
        love.graphics.setColor(INK)
        love.graphics.rectangle("line", gelX - lw - 3.5, my - 8.5, lw, 18)
        text(label, gelX - lw, my - 7, fonts.monoSmall, INK)
    end

    if not lanes[focus].ladder then drawAnnotations(lanes[focus], gelX + gelW, top, bottom) end

    local legendX = gelX + gelW + 100
    drawLegend(legendX, gelTop - 22, x + w - legendX)
end

local function niceStep(length, count)
    local raw = length / count
    local power = 10 ^ math.floor(math.log10(raw))
    for _, m in ipairs({ 1, 2, 5, 10 }) do
        if m * power >= raw then return m * power end
    end
end

local function drawMap(x, y, w, h)
    heading("D", "Linear restriction map", x, y, w, "hover a site to see its position")
    local tracks = {}
    for _, enzyme in ipairs(ENZYMES) do
        if enzyme.enabled then tracks[#tracks + 1] = enzyme end
    end
    tracks[#tracks + 1] = combined
    if #tracks == 0 then
        text("No enzyme selected.", x, y + h / 2 - 8, fonts.caption, GREY, "center", w)
        return
    end

    local left, right, length = x + 96, x + w - 64, #genome.dna
    local rulerY = y + 52
    local step = niceStep(length, 8)
    text("Position (bp)", x, rulerY - 8, fonts.small, GREY)
    love.graphics.setColor(INK)
    love.graphics.setLineWidth(1)
    love.graphics.line(left, rulerY, right, rulerY)
    for value = 0, length, step / 5 do
        local tx = left + value / length * (right - left)
        local major = value % step == 0
        love.graphics.setColor(INK)
        love.graphics.line(tx, rulerY, tx, rulerY - (major and 6 or 3))
        if major then text(commas(value), tx - 40, rulerY - 22, fonts.monoSmall, INK, "center", 80) end
    end

    local top = rulerY + 12
    local rowHeight = math.min(30, (y + h - top) / #tracks)
    local mx, my = love.mouse.getPosition()
    for i, track in ipairs(tracks) do
        local cy = top + (i - 0.5) * rowHeight
        if track.color then
            love.graphics.setColor(track.color)
            love.graphics.rectangle("fill", x, cy - 4, 8, 8)
        end
        text(track.name, x + 14, cy - 8, fonts.body, INK)
        rule(left, cy, right - left, 1, GREY)
        local tick = rowHeight * 0.3
        for _, cut in ipairs(track.cuts) do
            local owner = track.owners and track.owners[cut] or track
            local cx = left + cut / length * (right - left)
            love.graphics.setColor(owner.color)
            love.graphics.rectangle("fill", cx - 1, cy - tick, 2, tick * 2)
            if math.abs(mx - cx) <= 3 and math.abs(my - cy) <= tick then
                hover = { string.format("%s  cut after %s", owner.name, commas(cut)), mx, my }
            end
        end
        text(plural(#track.cuts, "site"), right + 10, cy - 7, fonts.small, GREY)
    end
end

local function drawTooltip()
    if not hover then return end
    local label, mx, my = hover[1], hover[2], hover[3]
    local tw = fonts.monoSmall:getWidth(label) + 14
    local tx = math.min(mx + 12, love.graphics.getWidth() - tw - 8)
    love.graphics.setColor(1, 1, 1)
    love.graphics.rectangle("fill", tx, my + 14, tw, 22)
    love.graphics.setColor(INK)
    love.graphics.setLineWidth(1)
    love.graphics.rectangle("line", tx + 0.5, my + 14.5, tw, 22)
    text(label, tx + 7, my + 18, fonts.monoSmall, INK)
end

function love.load(arg)
    love.graphics.setBackgroundColor(1, 1, 1)
    makeTextures()
    fonts.title = loadFont("LiberationSerif-Bold.ttf", 24)
    fonts.heading = loadFont("LiberationSerif-Bold.ttf", 16)
    fonts.caption = loadFont("LiberationSerif-Italic.ttf", 13)
    fonts.body = loadFont("LiberationSans-Regular.ttf", 13)
    fonts.bodyBold = loadFont("LiberationSans-Bold.ttf", 13)
    fonts.small = loadFont("LiberationSans-Regular.ttf", 11)
    fonts.smallBold = loadFont("LiberationSans-Bold.ttf", 11)
    fonts.mono = loadFont("LiberationMono-Regular.ttf", 13)
    fonts.monoSmall = loadFont("LiberationMono-Regular.ttf", 12)
    loadGenome(arg[1] or love.filesystem.getSource() .. DEFAULT_GENOME)
end

function love.update(dt)
    progress = math.min(1, progress + dt / 0.9)
end

function love.draw()
    local w, h = love.graphics.getDimensions()
    targets, hover = {}, nil
    drawHeader(w)
    if not genome then return end

    local top, leftW = 92, 430
    local rightX = 24 + leftW + 32
    local rightW, gelH = w - rightX - 24, math.floor((h - top - 20) * 0.6)
    local enzymesBottom = drawEnzymes(24, top, leftW)
    drawFragments(24, enzymesBottom + 22, leftW, h - enzymesBottom - 42)

    love.graphics.setColor(RULE)
    love.graphics.setLineWidth(1)
    love.graphics.line(rightX - 16, top, rightX - 16, h - 20)
    love.graphics.line(rightX, top + gelH + 10, w - 24, top + gelH + 10)
    drawGel(rightX, top, rightW, gelH)
    drawMap(rightX, top + gelH + 24, rightW, h - top - gelH - 44)
    drawTooltip()
end

function love.mousepressed(mx, my, button)
    if button ~= 1 then return end
    for _, target in ipairs(targets) do
        if mx >= target.x and mx <= target.x + target.w and my >= target.y and my <= target.y + target.h then
            return target.action()
        end
    end
end

function love.wheelmoved(_, dy)
    local mx, my = love.mouse.getPosition()
    if list and mx >= list.x and mx <= list.x + list.w and my >= list.y and my <= list.y + list.h then
        scroll = math.max(0, math.min(list.max, scroll - dy * ROW * 3))
    end
end

function love.filedropped(file)
    loadGenome(file:getFilename())
end

function love.keypressed(key)
    if key == "escape" then love.event.quit() end
end

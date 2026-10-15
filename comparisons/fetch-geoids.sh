#!/bin/sh
# Downloads GeographicLib's geoid grids and its 500 000-point test set into data/ (about 230 MB,
# 600 MB unpacked), for the geoid comparison. It also writes data/egm96-5.png, the same egm96-5
# grid as a PNG, which is what the egm96 crate builds from (see .cargo/config.toml).
set -eu
cd "$(dirname "$0")"
mkdir -p data
cd data
base=https://downloads.sourceforge.net/project/geographiclib
for g in egm84-15 egm96-15 egm96-5 egm2008-5 egm2008-2_5 egm2008-1; do
    [ -f "geoids/$g.pgm" ] || curl -fsSL "$base/geoids-distrib/$g.tar.bz2" | tar xjf -
done
[ -f GeoidHeights.dat ] || curl -fsSL "$base/testdata/GeoidHeights.dat.gz" | gunzip > GeoidHeights.dat
[ -f egm96-5.png ] || python3 - <<'EOF'
# A 16-bit greyscale PNG with exactly the PGM's pixels (both store them big-endian).
import struct, zlib
data = open("geoids/egm96-5.pgm", "rb").read()
fields, pos = [], 0
while len(fields) < 4:  # P5, width, height, maxval, skipping # comments
    while data[pos:pos + 1].isspace():
        pos += 1
    if data[pos:pos + 1] == b"#":
        pos = data.index(b"\n", pos)
        continue
    end = pos
    while not data[end:end + 1].isspace():
        end += 1
    fields.append(data[pos:end])
    pos = end
w, h = int(fields[1]), int(fields[2])
pixels = data[pos + 1:]
assert len(pixels) == 2 * w * h
rows = b"".join(b"\0" + pixels[2 * w * r:2 * w * (r + 1)] for r in range(h))
chunk = lambda kind, body: struct.pack(">I", len(body)) + kind + body + struct.pack(">I", zlib.crc32(kind + body))
png = (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 16, 0, 0, 0, 0))
       + chunk(b"IDAT", zlib.compress(rows, 9)) + chunk(b"IEND", b""))
open("egm96-5.png", "wb").write(png)
EOF
echo "geoid data ready in $(pwd)"

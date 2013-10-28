#!/usr/bin/env python
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
import os
import sys
import subprocess
from collections import defaultdict

working_dir = sys.argv[1]

for filename in os.listdir(sys.argv[1]):
    if filename.endswith(".png"):
        os.remove(os.path.join(working_dir, filename))


for filename in os.listdir(sys.argv[1]):
    if not filename.endswith(".dot"):
        continue

    fullname = os.path.join(working_dir, filename)
    subprocess.call("dot -Tpng " + fullname + " > " + fullname + ".png", shell=True)


compdict = defaultdict(list)

for filename in os.listdir(sys.argv[1]):
    if not filename.endswith(".png"):
        continue

    comp = filename.split("-")[0]
    compdict[comp].append(filename)

for component, files in compdict.iteritems():
    bg_name = component + "-bg.dot.png"
    prelinks_name = component + "-prelinks.dot.png"
    trees_names = filter(lambda x: x not in [bg_name, prelinks_name], files)
    trees_names = map(lambda t: os.path.join(working_dir, t), trees_names)

    bg = Image.open(os.path.join(working_dir, bg_name))
    prelinks = Image.open(os.path.join(working_dir, prelinks_name))
    trees = map(Image.open, trees_names)
    for t in trees:
        t.thumbnail((t.size[0] * 2 / 3, t.size[1] * 2 / 3), Image.ANTIALIAS)

    offset = max(bg.size[0], prelinks.size[0])
    width = sum(map(lambda s: s.size[0], trees)) + offset
    tree_height = max(map(lambda s: s.size[1], trees))
    height = max(bg.size[1] + prelinks.size[1], tree_height)

    #font = ImageFont.truetype("sans.ttf", 16)

    out = Image.new("RGB", (width, height), (255, 255, 255))
    out.paste(bg, (0, 0))
    out.paste(prelinks, (0, bg.size[1]))

    for tree, name in zip(trees, trees_names):
        draw = ImageDraw.Draw(out)
        name = os.path.basename(name)
        draw.text((offset, tree.size[1]), name, (0, 0, 0))
        out.paste(tree, (offset, 0))
        offset += tree.size[0]

    out.save(os.path.join(working_dir, component) + ".png")

    os.remove(os.path.join(working_dir, bg_name))
    os.remove(os.path.join(working_dir, prelinks_name))
    for t in trees_names:
        os.remove(t)

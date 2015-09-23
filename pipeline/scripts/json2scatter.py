#! /usr/bin/env python
#@Author Jose Fernandez
""" Script for creating a quality scatter plot from a json ST-data file.
The output will be a .png file with the same name as the json file stored in
the current directory.

If a regular expression for a gene symbol to highlight is provided, the output
image will be stored in a file called *.0.png. It is possible to give several
regular expressions to highlight, by adding another --highlight parameter.
The images from these queries will be stored in files with the number
increasing: *.0.png, *.1.png, *.2.png, etc.
"""

import argparse
import os
import re
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt
from stpipeline.common.json_utils import json_iterator
from stpipeline.common.json_utils import write_json
import numpy as np

blues = LinearSegmentedColormap.from_list('blues', [[0.02, 0.44, 0.69, 0.3],
                                                    [0.02, 0.44, 0.69, 1.0]])

def main(json_file, highlight_regexes, image, only_highlight, cutoff):
    fig = []
    ax = []
    highlights = []
    for _ in highlight_regexes if highlight_regexes else [0]:
        f = plt.figure()
        a = f.add_subplot(111, aspect='equal')
        fig.append(f)
        ax.append(a)
        highlights.append(set([]))

    expression = np.zeros((1000, 1000), dtype=np.int)

    it = json_iterator(json_file)
    for doc in it:
        
        if cutoff is not None and int(doc['hits']) < cutoff:
            continue
        
        expression[doc['x'], doc['y']] += doc['hits']

        if not highlight_regexes:
            continue

        for i, regex in enumerate(highlight_regexes):
            if re.search(regex, doc["gene"]):
                highlights[i].add((doc['x'], doc['y']))

    x, y = expression.nonzero()

    cmap = blues
    if image:
        img = np.load(image)
        img_edge = np.nonzero(img)[0].min()
        cmap = cm.Blues_r

    for a in ax:
        if not only_highlight:
            a.scatter(x, y, c=expression[x, y],
                      edgecolor="none",
                      s=10,
                      label="Expression",
                      cmap=cmap,
                      )

        if image:
            a.imshow(img)

    for i, regex in enumerate(highlight_regexes if highlight_regexes else [0]):
        if len(highlights[i]) == 0:
            continue

        x, y = zip(*highlights[i])
        ax[i].scatter(x, y, c="#CA0020",
                      edgecolor="#CA0020",
                      s=20,
                      label=highlight_regexes[i])

    for i, a in enumerate(ax):
        if not image:
            a.invert_yaxis()
        else:
            a.set_ybound(lower=img_edge)

        a.set_xlabel("X")
        a.set_ylabel("Y")
        a.legend()
        a.set_title("Scatter", size=20)

        if highlight_regexes:
            ending = ".{0}.png".format(i)
        else:
            ending = ".png"

        img_file = os.path.basename(json_file).replace(".json", ending)
        if image:
            fig[i].set_size_inches(16, 16)
        else:
            fig[i].set_size_inches(10, 8)

        fig[i].savefig(img_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("json_file", help="JSON ST-data file")
    parser.add_argument("--highlight", help="Regular expression for \
                        gene symbols to highlight in the quality \
                        scatter plot. Can be given several times.",
                        default=None,
                        type=str,
                        action='append')
    parser.add_argument("--image", default=None)
    parser.add_argument("--only-highlight", default=False, action='store_true')
    parser.add_argument("--cutoff", type=int, default=None)
    
    args = parser.parse_args()

    main(args.json_file, args.highlight, args.image, args.only_highlight, args.cutoff)

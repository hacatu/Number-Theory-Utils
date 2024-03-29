#!/usr/bin/env python3
import sys
from math import floor, ceil

if len(sys.argv) == 1:
	#print("https://img.shields.io/badge/coverage-n%2Fa-inactive")
	print("""<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="90" height="20" role="img" aria-label="coverage: n/a"><title>coverage: n/a</title><linearGradient id="s" x2="0" y2="100%"><stop offset="0" stop-color="#bbb" stop-opacity=".1"/><stop offset="1" stop-opacity=".1"/></linearGradient><clipPath id="r"><rect width="90" height="20" rx="3" fill="#fff"/></clipPath><g clip-path="url(#r)"><rect width="61" height="20" fill="#555"/><rect x="61" width="29" height="20" fill="#9f9f9f"/><rect width="90" height="20" fill="url(#s)"/></g><g fill="#fff" text-anchor="middle" font-family="Verdana,Geneva,DejaVu Sans,sans-serif" text-rendering="geometricPrecision" font-size="110"><text aria-hidden="true" x="315" y="150" fill="#010101" fill-opacity=".3" transform="scale(.1)" textLength="510">coverage</text><text x="315" y="140" transform="scale(.1)" fill="#fff" textLength="510">coverage</text><text aria-hidden="true" x="745" y="150" fill="#010101" fill-opacity=".3" transform="scale(.1)" textLength="190">n/a</text><text x="745" y="140" transform="scale(.1)" fill="#fff" textLength="190">n/a</text></g></svg>""")
else:
	perc_cov = float(sys.argv[1])
	red = floor((1 - perc_cov/100)*255)
	green = ceil(perc_cov/100*255)
	#print("https://img.shields.io/badge/coverage-{}%25-{:x}{:x}00".format(perc_cov, red, green))
	print(f"""<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="108" height="20" role="img" aria-label="coverage: {perc_cov}%"><title>coverage: {perc_cov}%</title><linearGradient id="s" x2="0" y2="100%"><stop offset="0" stop-color="#bbb" stop-opacity=".1"/><stop offset="1" stop-opacity=".1"/></linearGradient><clipPath id="r"><rect width="108" height="20" rx="3" fill="#fff"/></clipPath><g clip-path="url(#r)"><rect width="61" height="20" fill="#555"/><rect x="61" width="47" height="20" fill="#{red:x}{green:x}00"/><rect width="108" height="20" fill="url(#s)"/></g><g fill="#fff" text-anchor="middle" font-family="Verdana,Geneva,DejaVu Sans,sans-serif" text-rendering="geometricPrecision" font-size="110"><text aria-hidden="true" x="315" y="150" fill="#010101" fill-opacity=".3" transform="scale(.1)" textLength="510">coverage</text><text x="315" y="140" transform="scale(.1)" fill="#fff" textLength="510">coverage</text><text aria-hidden="true" x="835" y="150" fill="#010101" fill-opacity=".3" transform="scale(.1)" textLength="370">{perc_cov}%</text><text x="835" y="140" transform="scale(.1)" fill="#fff" textLength="370">{perc_cov}%</text></g></svg>""")


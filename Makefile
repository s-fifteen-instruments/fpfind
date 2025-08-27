PACKAGE=fpfind
DIR=src/${PACKAGE}
LIB=${DIR}/lib

.SILENT: shell
.PHONY: clean,read,readx

all:
	gcc -o ${DIR}/freqcd ${DIR}/freqcd.c ${LIB}/getopt.c -lm
shell:
	-poetry shell
fpfind-shell:
	-poetry run python -ic "import fpfind"

# Package-related stuff
test:
	poetry run pytest
test64:
	poetry run pytest --float64
clean:
	find . -type d -name ".pytest_cache" | xargs rm -rf
	find . -type d -name "__pycache__" | xargs rm -rf
	rm -f ${DIR}/freqcd
compile-schematic:
	java -jar docs/plantuml* docs/schematic.wsg

# Python package management
install: install-poetry
	poetry install
install-poetry:
	pip install -U poetry
uninstall:
	poetry run pip uninstall ${PACKAGE}

# Other shortcuts:
#
# - Convert 'a1' to 'a2':
#     xxd -e -g 4 -c 8 <FILE> | cut -b11-18,20-27
#
# - Convert 'a1' to 'a0':
#     xxd -e -g 4 -c 4 <FILE> | cut -b11-18 | sed -n 'h;n;p;g;p'

# ---------------------------------------------------------------------
# OS parsing

ifeq ($(OS),Windows_NT)
	OSFLAGS = -shared
	GCC = x86_64-w64-mingw32-g++.exe
	OUT = src/build/manyiv_windows.plugin
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		OSFLAGS = -shared -fPIC -DSYSTEM=OPUNIX
		OUT = src/build/manyiv_unix.plugin
	endif
	ifeq ($(UNAME_S),Darwin)
		OSFLAGS = -bundle -DSYSTEM=APPLEMAC
		OUT = src/build/manyiv_macosx.plugin
	endif
	GCC = g++
endif

ifeq ($(EXECUTION),windows)
	OSFLAGS = -shared
	GCC = x86_64-w64-mingw32-g++
	OUT = src/build/manyiv_windows.plugin
endif

EIGEN = /usr/include/eigen3
CFLAGS = -I $(EIGEN) -Wall -O3 $(OSFLAGS)

# ---------------------------------------------------------------------
# Main

## Compile directory
all: clean manyiv

## Download latest OSX plugin
osx_plugins:
	wget https://raw.githubusercontent.com/mcaceresb/stata-manyiv/osx/manyiv_macosx.plugin
	mv -f manyiv_macosx.plugin src/build/manyiv_macosx.plugin

# ---------------------------------------------------------------------
# Rules

## Compile manyiv plugin
manyiv: src/plugin/manyiv.cpp src/plugin/stplugin.c
	$(GCC) $(CFLAGS) -o $(OUT) $^

.PHONY: clean
clean:
	rm -f $(OUT)

#######################################################################
#                                                                     #
#                    Self-Documenting Foo (Ignore)                    #
#                                                                     #
#######################################################################

.DEFAULT_GOAL := show-help

.PHONY: show-help
show-help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}' \
	| more $(shell test $(shell uname) = Darwin && echo '--no-init --raw-control-chars')


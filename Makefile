CC		= gcc
CFLAGS		= -Wall -std=c18 -pedantic -Werror -Ofast
LDFLAGS		=
LDLIBS		= -lm -lgsl
ARCH		:= $(shell uname -m)
ifeq ($(ARCH),x86_64)
	LDLIBS	+= -lmvec -lcblas
endif

MAIN_SRC_DIR	= src/fit_executable
CONVERT_SRC_DIR	= src/convert_executable
MAIN_MODELS_DIR	= src/fit_executable/models
BINDIR		= bin

MAIN_SRC	= $(MAIN_SRC_DIR)/get_args.c \
		  $(MAIN_SRC_DIR)/get_r_squared.c \
		  $(MAIN_SRC_DIR)/read_file.c \
		  $(MAIN_SRC_DIR)/fit.c \
		  $(MAIN_SRC_DIR)/check.c \
		  $(MAIN_SRC_DIR)/analyze.c \
		  $(MAIN_SRC_DIR)/main.c \
		  $(MAIN_MODELS_DIR)/zdanovskii.c \
		  $(MAIN_MODELS_DIR)/uniquac.c \
		  $(MAIN_MODELS_DIR)/norrish.c \
		  $(MAIN_MODELS_DIR)/virial.c \
		  $(MAIN_MODELS_DIR)/caurie.c \
		  $(MAIN_MODELS_DIR)/raoult.c

CONVERT_SRC	= $(CONVERT_SRC_DIR)/get_args.c \
		  $(CONVERT_SRC_DIR)/read_file.c \
		  $(CONVERT_SRC_DIR)/write_to_file.c \
		  $(CONVERT_SRC_DIR)/convert.c \
		  $(CONVERT_SRC_DIR)/main.c \

MAIN_DEPS	= $(MAIN_SRC_DIR)/definitions_and_headers.h
MAIN_EXECUTABLE	= $(BINDIR)/FitWaterActivity
MAIN_SYS_FILE	= /usr/bin/FitWaterActivity

CONVERT_DEPS		= $(CONVERT_SRC_DIR)/definitions_and_headers.h
CONVERT_EXECUTABLE	= $(BINDIR)/ConvertWaterActivity
CONVERT_SYS_FILE	= /usr/bin/ConvertWaterActivity

.PHONY: all clean install uninstall
all: $(MAIN_EXECUTABLE) $(CONVERT_EXECUTABLE)

clean:
	rm -f $(MAIN_EXECUTABLE) $(CONVERT_EXECUTABLE)
	rmdir bin

install: $(MAIN_EXECUTABLE) $(CONVERT_EXECUTABLE)
	cp $(MAIN_EXECUTABLE) $(MAIN_SYS_FILE)
	cp $(CONVERT_EXECUTABLE) $(CONVERT_SYS_FILE)

uninstall: $(MAIN_SYS_FILE) $(CONVERT_SYS_FILE)
	rm -f $(MAIN_SYS_FILE) $(CONVERT_SYS_FILE)

$(MAIN_EXECUTABLE): $(MAIN_SRC) $(MAIN_DEPS)
	if [ ! -d "./bin" ]; then mkdir bin; fi
	$(CC) $(CFLAGS) $(LDFLAGS) $(LDLIBS) $(MAIN_SRC) -o $@

$(CONVERT_EXECUTABLE): $(CONVERT_SRC) $(CONVERT_DEPS)
	if [ ! -d "./bin" ]; then mkdir bin; fi
	$(CC) $(CFLAGS) $(LDFLAGS) $(LDLIBS) $(CONVERT_SRC) -o $@



TARGET = la-tests

CC = gcc
INCLUDE_DIR = src/LA

CFLAGS = -I$(INCLUDE_DIR) -I$(INCLUDE_DIR)/tests -O3
C_SOURCES = $(wildcard $(INCLUDE_DIR)/*.c $(INCLUDE_DIR)/tests/*.c)

DEPS = $(wildcard $(INCLUDE_DIR)/*.h $(INCLUDE_DIR)/tests/*.h)
OBJ = ${C_SOURCES:.c=.o}

# First rule is the one executed when no parameters are fed to the Makefile


%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(TARGET): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

run: $(TARGET)
	./$(TARGET)

clean:
	$(RM) LA/*.bin LA/*.o LA/*.dis LA/*.elf
	$(RM) lib/*.o

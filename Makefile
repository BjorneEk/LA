TARGET = LA
TEST_TARGET = LA_test
CC = gcc
INCLUDE_DIR = src
CFLAGS = -I$(INCLUDE_DIR) -F /Library/Frameworks -F /opt/homebrew/opt/openblas/lib -I /opt/homebrew/opt/openblas/include
AR_FLAGS = ruv
LIBS = -framework UL  -L/opt/homebrew/opt/openblas/lib -lopenblas -lpthread
C_TEST_SOURCES = $(wildcard src/*.c *.c src/*/*.c src/tests/test/*.c )
C_SOURCES = $(wildcard src/*/*.c src/*.c *.c )
DEPS = $(wildcard $(INCLUDE_DIR)/*/*.h $(INCLUDE_DIR)/*.h *.h)

TEST_OBJ = ${C_TEST_SOURCES:.c=.o}
OBJ = ${C_SOURCES:.c=.o}


FRAMEWORK_DIR = /Library/Frameworks/LA.framework
HEADER_DIR = /usr/local/include/LA

# First rule is the one executed when no parameters are fed to the Makefile


%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

$(TEST_TARGET): $(TEST_OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

run: $(TEST_TARGET)
	./$(TEST_TARGET)

framework: $(OBJ)
	ar $(AR_FLAGS) $(TARGET) $^

# updates the library framework in the framework dir and include dir
update: framework
	sudo cp $(DEPS) $(FRAMEWORK_DIR)/Headers
	sudo cp $(DEPS) $(HEADER_DIR)
	sudo cp $(TARGET) $(FRAMEWORK_DIR)

clean:
	$(RM) src/*.bin src/*.o src/*.dis src/*.elf
	$(RM) lib/*.o
	$(RM) $(TARGET)

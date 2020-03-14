PKG_CONFIG_PKGS = sdl glib-2.0
CFLAGS  = -W -Wall -O3 -std=gnu99 -march=native -mtune=native $(shell pkg-config --cflags $(PKG_CONFIG_PKGS))
LDFLAGS = 
LDLIBS  = -lm $(shell pkg-config --libs $(PKG_CONFIG_PKGS))

TARGET = smoke
OBJS = $(TARGET).o

.PHONY: all clean

all: $(TARGET)

clean:
	$(RM) $(TARGET) *.o

$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)


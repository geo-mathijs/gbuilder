MODULE_big = gbuilder
PG_CPPFLAGS = --std=c99
OBJS = $(patsubst %.c,%.o, $(wildcard src/*.c))
SHLIB_LINK += -llwgeom -lpredicates
EXTENSION = gbuilder
DATA = $(wildcard sql/*.sql)

PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)

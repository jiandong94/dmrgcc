# all subdirs
SUBDIRS=$(shell ls -l | grep ^d | awk '{print $$9}')

# .cc and .o files at current path
CUR_DIR=$(shell pwd)
CUR_SOURCE=$(wildcard *.cc)
CUR_OBJ=$(patsubst %.cc, $(OBJ_DIR)/%.o,$(notdir $(CUR_SOURCE)))

all:$(SUBDIRS) $(CUR_OBJ)

$(SUBDIRS):ECHO
	make -C $@
$(CUR_OBJ):$(OBJ_DIR)/%.o:$(CUR_DIR)/%.cc
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBRARY) -c $< -o $@
ECHO:
	@echo $(SUBDIRS)


# .cc and .o files at current path
CUR_SOURCE=$(wildcard *.cc)
#CUR_OBJ=$(patsubst %.cc, $(OBJ_DIR)/%.o, $(notdir $(CUR_SOURCE)))
#ALL_OBJ=$(wildcard $(OBJ_DIR)/*.o)
TARGET=$(patsubst %.cc, $(BIN_DIR)/%, $(notdir $(CUR_SOURCE)))
#CUR_BIN = test_real_matrix
all:$(TARGET) 

$(TARGET):$(BIN_DIR)/%:%.cc
	$(CC) $< $(OBJ_DIR)/*.o $(CFLAGS) $(INCLUDE) $(LIBRARY) -o $@
#$(CUR_OBJ):$(OBJ_DIR)/%.o:%.cc
#	$(CC) $(CFLAGS) $(INCLUDE) $(LIBRARY) -c $< -o $@

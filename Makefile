SRCS = inter_work_decomp.c work_decomp_tool.c tool_func.c
OUT_DO = do_dw_work_decomp_ala

$(OUT_DO): $(SRCS:.c=.o)
	gcc $(SRCS:.c=.o) -o $(OUT_DO) -lm

$(SRCS:.c=.o): $(SRCS)
	gcc -c $(SRCS)

#define TICK_AND_PRINT(counter, period, str) do { \
    if ((counter) == 0) { \
        hw.PrintLine(str); \
        (counter) = 1; \
    } else { \
        (counter)++; \
        if ((counter) == (period)) { \
            (counter) = 0; \
        } \
    } \
} while (0)


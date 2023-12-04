data = {
            0    : 1,
            0.25 : 0.989616,
            0.5  : 0.958851,
            0.75 : 0.908852,
            1    : 0.841471,
            1.25 : 0.759188,
            1.5  : 0.664997,
            1.75 : 0.562278,
            2    : 0.454649
        }


def TrapezoidMethod(data : dict[float, float]):
    x = list(data.keys())
    step = x[1] - x[0]
    values = list(data.values())
    sum = (values[0] + values[len(values) - 1]) / 2
    
    for i in range(1, len(values) - 1):
        sum += values[i]

    return sum * step


def SimpsonMethod(data : dict[float, float]):
    x = list(data.keys())
    step = x[1] - x[0]
    values = list(data.values())
    sum = values[0] + values[len(values) - 1]

    for i in range(1, len(values) - 1):
        sum += (2 + 2 * (i % 2 == 1)) * values[i] 

    return sum * step * 1 / 3


def main():
    items = list(data.items())
    data_2h = {}
    for i in range(0, len(data), 2):
        item = items[i]
        data_2h[item[0]] = item[1]

    I_h = TrapezoidMethod(data)
    I_2h = TrapezoidMethod(data_2h)
    
    I_sim = SimpsonMethod(data)
    I_r = I_h + (I_h - I_2h) / (2**2 - 1)

    print(f'I_h = {I_h}')
    print(f'I_2h = {I_2h}')
    print(f'I_sim = {I_sim}')
    print(f'I_r = {I_r}')


if __name__ == '__main__':
    main();

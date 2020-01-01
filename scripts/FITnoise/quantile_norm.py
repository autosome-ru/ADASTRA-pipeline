import numpy as np


def between(a, b, ratio):
    return (1 - ratio) * a + ratio * b


def interpolate_(values, reference):
    if len(reference) == 1:
        return [reference[0]] * len(values)
    lenv = values[-1] - values[0]
    segments = len(reference) - 1
    interpolated_values = []
    for value in values:
        if lenv == 0:
            part = 0.5
        else:
            part = (value - values[0]) / lenv
        index = int(part * segments)
        a = reference[index]
        b = reference[min(a + 1, len(reference) - 1)]
        ratio = (part - index / segments) * segments
        interpolated_values.append(between(a, b, ratio))
    return interpolated_values


def interpolate(values, reference):
    if len(reference) == 1:
        return [reference[0]] * len(values)
    lenv = len(values) - 1
    segments = len(reference) - 1
    interpolated_values = []
    for i in range(len(values)):
        if lenv == 0:
            part = 0.5
        else:
            part = i / lenv
        index = int(part * segments)
        a = reference[index]
        b = reference[min(a + 1, len(reference) - 1)]
        ratio = (part - index / segments) * segments
        interpolated_values.append(between(a, b, ratio))
    return interpolated_values


def rank_interpolation(values, reference_lists, coeffs):
    interpolated_values = np.zeros(len(values), dtype=np.float64)
    for reference, coef in zip(reference_lists, coeffs):
        print(coef, interpolate(values, reference))
        interpolated_values += coef * np.array(interpolate(values, reference))
    return interpolated_values


def get_lists_and_coefs(labeled_values, labeled_reference, unique_values):
    lists_and_coefs = []
    for value_count in unique_values:
        print(value_count)
        lists = []
        coefs = []
        inds, values = zip(
            *[(index, value) for index, (value, count) in enumerate(labeled_values) if count == value_count])
        amount = len(values)
        corresponding_reference_counts = [labeled_reference[index][1] for index in inds]
        for ref_count in list(set(corresponding_reference_counts)):
            coefs.append(sum(1 for count in corresponding_reference_counts if count == ref_count) / amount)
            lists.append([value for value, count in labeled_reference if count == ref_count])
        lists_and_coefs.append((values, lists, coefs))
    return lists_and_coefs


def normalize_values(full_values, full_reference, values_counts, reference_counts):
    order = np.argsort(np.array(list(zip(full_values, values_counts)), dtype=[('x', np.int_), ('y', np.int_)]),
                       order=('y', 'x')).argsort()
    print(order)
    labeled_values = sorted(list(zip(full_values, values_counts)), key=lambda x: x[0])
    labeled_values = sorted(list(zip(full_values, values_counts)), key=lambda x: x[1])
    labeled_reference = sorted(list(zip(full_reference, reference_counts)), key=lambda x: x[0])
    labeled_reference = sorted(list(zip(full_reference, reference_counts)), key=lambda x: x[1])
    unique_value_counts = sorted(list(set(values_counts)))
    normalized_values = []
    for values, lists, coefs in get_lists_and_coefs(labeled_values, labeled_reference, unique_value_counts):
        print(values, lists, coefs)
        print(rank_interpolation(values, lists, coefs))
        normalized_values.extend(rank_interpolation(values, lists, coefs))
    return np.round(normalized_values).astype(int)[order]


if __name__ == '__main__':
    vals = [10, 11, 12, 13, 15, 20, 25, 27, 38, 42]
    normed = np.array(
        normalize_values(vals, vals, [40, 30, 19, 9, 0, 2, 0, 1, 1, 1], [42, 31, 12, 5, 1, 0, 1, 1, 0, 0]))
    normed = np.round(normed).astype(int)
    print(normed)

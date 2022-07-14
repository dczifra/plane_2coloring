# Finding maximal 1-avoiding distance on plane
The goal of this repository is to find as large 1-avoiding set in the plane, as possible (largest density ==> *d*).

A 1-avoiding set is a set in the plane, where neither of the points in the set have 1-distance on plane.

Examples:
* Empty set: ==> *d({Empty}) = 0*
* Croft's construction ==> *d({Croft's construction}) = 0.22936*

# TODO
* square grid
* output indexes to file
* test create graph

# Result

## Square grid
```
|   R	|  Size	|Density (%)| Time  |
|-------|-------|-------    |-------|
|  0.02 |800x800|  19.7     | 38s   |
|  0.01 |800x800|  18.9     | 60s   |
|  0.005|800x800|  18.9     | 100s  |

|  0.02 |1.6x1.6|  19.6     | 2min  |
|  0.01 |1.6x1.6|  20.1     | 5min  |
|  0.005|1.6x1.6|  19.2     | 5min  |

|  0.01 |3.2x3.2|  20.2     | 10min  |
|  0.005|3.2x3.2|   -       | 10min  |
| 0.0025|3.2x3.2|   -       | 10min  |
```

## Visualization

The 800x800 grid:

![Torus with 800x800 grid ==> 19.7%"](/doc/grid.png)

The Hexagon grid:

![Hexagon grid based torus radius=200](/doc/hex1.png)



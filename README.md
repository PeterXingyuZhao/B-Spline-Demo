# B-Spline Demo

This project is about implementing the B-Spline curves and demonstrate their key features intuitively, as well as providing a concise introduction to the B-Spline concept.

In this project, I calculated the blending functions symbolically using the Cox-de Boor recursion algorithm and stored all the numbers as fractions for later display purpose.

The key features I demonstrated here are the degree, control points, knot vector, and the blending functions. I provided an interface where the user can set the degree of the polynomial blending functions, add/delete/move the control points, and specify the knot vector they want. (The default is the open uniform knot vector.) Based on these inputs, my program will calculate the blending functions, and display both the piecewise mathematical representation and their corresponding graphs.

In addition, the demo also provides a feature that the user can highlight a local region associated with a specific u value. It will highlight the control points and the blending functions which are effective.

The code uses the framework (including the index page structure and draggable points), which was primarily developed by Prof. Michael Gleicher with assistance from the CS 559 course staff, including Young Wu, over the years.

The content is Copyright &copy; 2025, Peter Xingyu Zhao.

This workbook is provided under a Creative Commons Attribution-NonCommercial 4.0 International license. See https://creativecommons.org/licenses/by-nc/4.0/ for the explanation and https://creativecommons.org/licenses/by-nc/4.0/legalcode for the license itself.

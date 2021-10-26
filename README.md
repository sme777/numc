# Welcome to NumC (like NumPy but by C)


## Matrix functions in C
###### allocate_matrix
``allocate_matrix`` was actually one of the last functions we wrote in this part. We thought the implementation would be slightly harder as opposed to other functions we needed to fill, so we left it to the end. Turned out it wasn't as hard as we thought after all.
###### allocate_matrix_ref
``allocate_matrix_ref`` was the last function we worte in this part 1. Also, it was the primary cause of pain and agony in the debugging. Initially, we thought that the differences between ``allocate_matrix`` and this function were the offsets and the exisiting matrix. We did not realize that the underlying structure need be preserved in both matricies. we fixed this in the debugging process of the Python-C interface, and although in the end it seemed quite straightforward, the pointer arithmetic was quite difficult at first. 
###### deallocate_matrix
We implemented ``deallocate_matrix`` function right after we wrote ``allocate_matrix``. It made sense to us to do it that way, cause once we understood the general structure of how the matrix is assigned and saved on the heap, we could easily deallocate it.
###### get
Writing the get function took a relatively short time, it was basically returning a specific indexed value.
###### set
Writing the set function also took relatively hsort time, since much like the get function, it was essentially referencing a nested array and setting a value.
###### fill_matrix
The code of the ``fill_matrix`` was a mixture of the add_matrix and set functions. We implemented this function by running 2 for loops, and essentially calling the set function with each index values.
###### add_matrix
Adding two matricies was the first function we implemented, we thought it was conceptionally easier to get a grasp on, and in the process we could understand how the approach the problem. It turned out to be an easy function to implement once the structure of the matrix was clear.
###### sub_matrix
This was the second function we implemented. Much of the code assembles that of the ``add_matrix``, with the only difference in the subtraction part.
###### mul_matrix
``mul_matrix`` turned out to take the most of the time in this part. Even with that much spent, we ended up getting it wrong once we completed Task 1. Albeit passing the test, we did not check whether ``mul_matrix`` worked on large matricies, and this was cause of mcuh pain later on. We ended up fixing ``mul_matrix`` in the debugging process of Python-C interface.
###### pow_matrix
Writing the ``pow_matrix`` function also took some time, although it wasn't so hard conceptually as the ``mul_matrix``, as we were calling the mul_matrix to handle the heavy lifting. One thing we missed out was the case when the power is taken of 0 or 1. These two cases were handled later on in the debugging process of the Python-C interface part.
###### neg_matrix
Negating a matrix seemed to be not too difficult of a task, so we did that third, once we was done with the ``add_matrix`` and ``sub_matrix``.
###### abs_matrix
Taking the absolute of the function was quite similar to the ``neg_matrix`` function, the only difference was to include the conditional to see if a value at matrix[i][j] was negative.

## Writing the Python-C interface

#### Number Methods
The number methods turned out to be the easiest of the all function to write, and also took the least time. Here is detailed description for each function.
###### Matrix61c_add
It took us some time to understand how the Python-C interface would work and the correct way to implement this function. Essentially, it was calling the function we already wrote in part 1, with the caveat that we needed to unpack arguments and do error checking. After reading the Python docs on the C interface, we implemented this function with ease.
###### Matrix61c_sub
Matrix61c_sub was much like ``Matrix61c_add``, we still needed to do type casting of python objects and also error checking, after which we just called the matrix_sub function we already implemented in part 1.
###### Matrix61c_multiply
While initially it seemed that ``Matrix61c_multiply`` would be very different from Matrix61c_add and others, it turned out to be very similar. We still had to do type checks, castings and error handling.
###### Matrix61c_neg
This function was the first one we implemented. Although it turned out not to be too different from ``Matrix61c_sub`` or ``Matrix61c_add``, we thought it was a good starting point since we needed to do less error checking and object casting.
###### Matrix61c_abs
The implementation for ``Matrix61c_abs`` was very similar to that of ``Matrix61c_neg``, with one difference in calling the matrix_abs funcion instead of ``neg_matrix``.
###### Matrix61c_pow
Much like the rest of the functions, we type checked python object before we casted and then did error handling both before and after calling ``pow_matrix`` function in matrix.c.
###### Matrix61c_as_number
At first, we were confused on how to deal with this functions, but after reading the docs we found the right way and got rid of compiler warnings.

#### Instance Methods
Instance methods took collectively less time to implement than the Number methods, but per function, we spent more time on the Instance methods. Below is a detailed description of each function.
###### Matrix61c_set_value
For Matrix61c_set_value, we first needed to check the number of argument supplied to the ``PyObject* args``, and make sure that it was 3. After, we basically did bunch of type checks and castings, eventually calling the ``set()`` function in the matrix.c. 
###### Matrix61c_get_value
Matrix61c_get_value was very similar to Matrix61c_set_value with the caveat that we had less values in the ``PyObject* args``. So, we made sure that the arguments supplied were 2, and did type checks, casting and finally called ``get()`` in matrix.cc 
###### Matrix61c_methods
Frankly, this was the most tedious part in the whole project. First because we had no idea what was needed to be done, and the docs were not really helpful. After finally finding a workaround for the Matrix61c_get_value function, the Matrix61c_set_value started failing and then we were very perplexed. Turns out in the number of arguments, we are not actually supposed to hardcode those values. Just solving this tiny issue took a lot of time.
#### Indexing
This was our sneak peek to hell. Writing the each of these functions and having them work properly took days.
###### Matrix61c_subscript
Writing ``Matrix61c_subscript`` took around 250 lines of code. Since there were a large amount of possible combinations passed to the ``PyObject* key``, each had to be checked seperately, both for 1D and 2D matricies. First, we implemented the subscript for 1D matricies, which took collectively about a day and a half. After thaving the 1D case work, we expanded to the general 2D case. The main difference between the 1D and 2D matricies for this marticular function weas having an extra type o variables for ``key``, in particular, the tuple of integer/slice. 
###### Matrix61c_set_subscript
Writing ``Matrix61c_set_subscript`` took considerably less time now that we had finished writing the ``Matrix61c_subscript`` function. Actually, writing ``Matrix61c_subscript`` first made it easier for us, since we can call this function in the implementation of ``Matrix61c_set_subscript``, and so that saved a lot of code lines. We still had the problem of configuring whether the passed arguments were for a 1D matrix. We resorted to making a helper function named ``isASingleNumberIndex`` which did just that. It took two argument ``Matrix61c *self`` and  ``PyObject *key`` and returned 1 if the matrix was 1D and 0 if not. After checking if a matrix is 1D we simply implemented the ``set()`` function in applicable places. 

As a side note, in both functions that fell under indexing we did an incredible amount of error and type checking. This made sure when the unreasonable arguments were passed, a proper error was spitted out.  

## Speedups

###### allocate_matrix
###### allocate_matrix_ref
###### deallocate_matrix
###### get
###### set
###### fill_matrix
###### add_matrix
###### sub_matrix
###### mul_matrix
###### pow_matrix
###### neg_matrix
###### abs_matrix
Making contiguous memory. This was done in the ``allocate_matrix`` and ``allocate_matrix_ref``. Then we used openmp on any matrices with more than 100k elements, this made sure that the smaller matricies were not affected by the openmp yet the larger matricies still yieled results pretty fast. Next step was to use repeated squaring for log(n) runtime in ``pow_matrix``. After, we unrolled for loops when appropriate, usually 4x, in function specified above. Finally, we used SIMD in ``mul_matrix``, and also loop unrolling and openmp, and transposed the mat2 for better cache hits.


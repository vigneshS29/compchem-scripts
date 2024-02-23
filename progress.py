#Author: Vignesh Sathyaseelan (vsathyas@purdue.edu)
import time,io

'''
Progress Bar Wrapper

 This Python script provides a flexible progress bar wrapper to enhance your iterative processes.

 Features

Adaptable: Works seamlessly with:
range objects for tracking progress over a known range of iterations.
General iterables (lists, custom iterators, etc.) for when the total length is unknown.
NumPy arrays to visualize progress during array processing.
Files to add progress indicators for file operations.
enumerate for iterations where both index and value are needed.

'''

def prange(iterable,total=None,static=None):
    start = time.time()
    if isinstance(iterable, enumerate):
        if total == None: raise ValueError ('Please Provide total in case of enumerate')
        bar_length = 50  # Adjust this for desired bar length
        for current_step,item in iterable:
            percent_complete = (current_step + 1) / total * 100
            filled_bar = int(bar_length * percent_complete / 100)
            progress_bar = "[" + "#" * filled_bar + "-" * (bar_length - filled_bar) + "]"
            if static: print(f" \r Progress (Enumerate) : {progress_bar} {percent_complete:.2f}% Time: {time.time() - start:.10f}s", end="") 
            else: print(f"Progress (Enumerate) : {progress_bar} {percent_complete:.2f}% Time: {time.time() - start:.10f}s \n", end="") 
            yield(current_step,item)

    elif isinstance(iterable, io.IOBase):
        iterable = iterable.readlines()
        total = len(iterable)
        bar_length = 50  # Adjust this for desired bar length
        for current_step,item in enumerate(iterable):
            percent_complete = (current_step + 1) / total * 100
            filled_bar = int(bar_length * percent_complete / 100)
            progress_bar = "[" + "#" * filled_bar + "-" * (bar_length - filled_bar) + "]"
            if static: print(print(f" \r Progress (File Read) : {progress_bar} {percent_complete:.2f}% Time: {time.time() - start:.10f}s", end="") )
            else:print(f"Progress (File Read) : {progress_bar} {percent_complete:.2f}% Time: {time.time() - start:.10f}s \n", end="") 
            yield(item)
      
    else:
        total = len(iterable)
        bar_length = 50  # Adjust this for desired bar length
        for current_step,item in enumerate(iterable):
            percent_complete = (current_step + 1) / total * 100
            filled_bar = int(bar_length * percent_complete / 100)
            progress_bar = "[" + "#" * filled_bar + "-" * (bar_length - filled_bar) + "]"
            if static: print(f" \r Progress (Iterate) : {progress_bar} {percent_complete:.2f}% Time: {time.time() - start:.10f}s", end="") 
            else: print(f"Progress (Iterate) : {progress_bar} {percent_complete:.2f}% Time: {time.time() - start:.10f}s \n", end="") 
            yield(item)

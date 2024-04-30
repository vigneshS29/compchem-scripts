#Author: Vignesh Sathyaseelan (vsathyas@purdue.edu)
import time,io

'''
Progress Bar Wrapper

This Python script provides a flexible progress bar wrapper to enhance your iterative processes.

Features:

Adaptable: Works seamlessly with:
range objects for tracking progress over a known range of iterations.
General iterables (lists, custom iterators, etc.) for when the total length is unknown.
NumPy arrays to visualize progress during array processing.
Files to add progress indicators for file operations.
enumerate for iterations where both index and value are needed.

'''

def io_progress(name, current_step, total, start, static, bar_length=50):
    
    percent_complete = (current_step + 1) / total * 100
    filled_bar = int(bar_length * percent_complete / 100)
    progress_bar = "[" + "#" * filled_bar + "-" * (bar_length - filled_bar) + "]"
    if static: print(f" \r Progress {name} : {progress_bar} {percent_complete:.2f}% Time: {time.time() - start:.10f}s", end="") 
    else: print(f"Progress {name} : {progress_bar} {percent_complete:.2f}% Time: {time.time() - start:.10f}s \n", end="")
    
    return 

def prange(iterable,total=None,static=None):
    
    start = time.time()
    
    if isinstance(iterable, enumerate):
        if total == None: raise ValueError ('Please Provide total in case of enumerate')
        for current_step,item in iterable:
            io_progress('Enumerate',current_step,total,start,static)
            yield(current_step,item)

    elif isinstance(iterable, io.IOBase):
        iterable = iterable.readlines()
        total = len(iterable)
        for current_step,item in enumerate(iterable):
            io_progress('File',current_step,total,start,static)
            yield(item)
      
    else:
        total = len(iterable)
        for current_step,item in enumerate(iterable):
            io_progress('Iterable',current_step,total,start,static)
            yield(item)
    return

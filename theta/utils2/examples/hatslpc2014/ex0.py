# IMPORTANT NOTE: python requires consistent indendation. Please configure
# your text editor to use spaces only as indentation, in units of 4 spaces per tab.

execfile("common.py")
print "Hello World!"

# example for "if" statement and indentation:
if True:
    print "Hello World from indented code."

    
# while loop example (will not execute):
while False:
    print "You will not see this."
    
# method definition example with default argument value:
def print_s(s = ''):
    # The '%' operator works similar to C's sprintf:
    print "called print_s with s=%s" % s


# method invocation will use default arguments:
print_s()
# method invocation with "argument name = value" notation
print_s(s = 'hello world')    

# make a list:
l = [1,2,3]

# access list items through 0-based indices:
print "list items 0, 2:", l[0], l[2]

# for loop over list items:
print "list items:"
for li in l:
    print li

    
# Test the theta setup and calculate a simple likelihood ratio 
model = test_model.simple_counting(s = 1.0, b = 10.0)
lrs = get_bkg_t(model)
plot_histogram(lrs, "test.pdf")

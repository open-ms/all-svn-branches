import sandbox as s

try:
    l = s.DoubleList([1,2,3, 4+2j])
except:
    import traceback
    traceback.print_exc()

dv = s.DataValue(map(float, [1,2,3]))
print dv.valueType()
try:
    dv = s.DataValue([1.0,2,1+2j])
except:
    import traceback
    traceback.print_exc()


mse = s.MSExperiment()
print mse.getMetaValue("hi")

fh = s.FileHandler()
#fh.loadExperiment("example1.mzML", mse)
fh.loadExperiment("example1.mzXML", mse)
print mse.getMetaValue("hi")
print mse.getStringKeys()
print mse.getIntKeys()
print mse.setMetaValue("hi", 4711)
print mse.getMetaValue("hi")
print mse.getStringKeys()
print mse.getIntKeys(["a"])



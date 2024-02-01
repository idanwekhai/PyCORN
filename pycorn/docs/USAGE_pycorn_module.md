# Using the pycorn module:

Import the Unicorn >= 6 interface:

```python
from pycorn import PcUni6
```

Create an instance of the pycorn object.

```python
my_res_file = PcUni6("sample1.res")
```

Parse the file. This will load all tabular data (e.g. UV, Cond, etc) that is stored as xml in the .res file:

```python
my_res_file.load_all_xml()
```

Alternatively: load all keys from the `.res` file and inspect them manually.

```python
my_res_file.load()
print(my_res_file.keys())

>>>[u'CreationNotes', u'Logbook', u'UV', u'Cond', u'pH', u'Pressure', u'Temp', u'Conc', u'Fractions']
```

The above list is your key to access the data inside the file as
`my_res_file[key][sub_key]`. ``sub_key`` can be:

- ``data``: contains the actual data, either pure text or a list with x/y-pairs as tuples
- ``unit``: the unit for this data block (mAu, ms/cm etc.)
- ``run_name``: an internal name (like "Manual Run 8")

For example to read-out the the unit for the UV-block:

```python
my_res_file['UV']['unit']
```

To read the actual UV-data:

```python
x = my_res_file['UV']['data']
print(x[0:3])
>>>[(0.0, -9.22), (0.06, -0.007), (0.13, -0.004)]
```

function result = mat_fcn(signal)
 model = py.importlib.import_module('FConvNetmin');
 py.importlib.reload(model);
 fcn = py.importlib.import_module('model_for_matlab');
 py.importlib.reload(fcn);
 result = fcn.process(pyargs('data', signal));
end
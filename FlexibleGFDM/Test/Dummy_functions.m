function func = Dummy_functions()
% Structure of reference to functions
func.func1 = @func1;
func.func2 = @func2;
% implementation of the functions
    function func1 ()
        disp('This is function 1');
    end
    function func2 ()
        disp('This is function 2');
    end
end

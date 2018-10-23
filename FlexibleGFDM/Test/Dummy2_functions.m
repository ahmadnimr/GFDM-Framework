function func = Dummy2_functions()
% Structure of reference to functions
func.func3 = @func3;
func.func4 = @func4;
% use functions from Dummy_functions
dum1 = Dummy_functions();
% implementation of the functions
    function func3 ()
        disp('This is a call of function Dummy func 1');
        dum1.func1();
    end
    function func4 ()
         disp('This is a call of function Dummy func 2');
         dum1.func2();
    end
end

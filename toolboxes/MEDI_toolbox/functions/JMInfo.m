classdef JMInfo
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        s
    end

    methods
        function obj = JMInfo(s)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.s=s;
        end

        function outputArg = nonzeros(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = [];
        end

        function outputArg = isreal(obj) %#ok<*MANU> 
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = true;
        end

        function out=isoptimargdbl(~, ~, ~, obj)
               out=[];
        end

        function s=sparse(obj)
            s=obj.s;
        end

    end
end
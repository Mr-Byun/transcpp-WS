%CLASS_INTERFACE Example MATLAB class wrapper to an underlying C++ class
classdef organism_interface < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance
        function this = organism_interface(varargin)
            this.objectHandle = organism_interface_mex('new', varargin{:});
        end
        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            organism_interface_mex('delete', this.objectHandle);
        end

        %% Test - example class method call
        function varargout = test(this, varargin)
            [A] = organism_interface_mex('test', this.objectHandle, varargin{:});
            [varargout{1:nargout}] = A;
        end
        
        %% getNNuc - get number of nuclei
        function varargout = getNNuc(this, varargin)
            n = organism_interface_mex('getNNuc', this.objectHandle, varargin{:});
            [varargout{1:nargout}] = double(n);
        end
        
        %% getNParams - get number of nuclei
        function varargout = getNParams(this, varargin)
            n = organism_interface_mex('getNParams', this.objectHandle, varargin{:});
            [varargout{1:nargout}] = double(n);
        end

        %% getNGenes - get number of genes
        function varargout = getNGenes(this, varargin)
            n = organism_interface_mex('getNGenes', this.objectHandle, varargin{:});
            [varargout{1:nargout}] = double(n);
        end

        %% getTotalScore - get model score
        function varargout = getTotalScore(this, varargin)
            [varargout{1:nargout}] = organism_interface_mex('getTotalScore', this.objectHandle, varargin{:});
        end

        %% getParameterName - get model score
        function varargout = getParameterName(this, varargin)
            [varargout{1:nargout}] = organism_interface_mex('getParameterName', this.objectHandle, varargin{:});
        end

        %% getParameterValue - get model score
        function varargout = getParameterValue(this, varargin)
            [varargout{1:nargout}] = organism_interface_mex('getParameterValue', this.objectHandle, varargin{:});
        end

        %% setParameterValue - get model score
        function varargout = setParameterValue(this, varargin)
            [varargout{1:nargout}] = organism_interface_mex('setParameterValue', this.objectHandle, varargin{:});
        end

        %% getN - get model Ns
        function varargout = getN(this, varargin)
            [varargout{1:nargout}] = organism_interface_mex('getN', this.objectHandle, varargin{:});
        end
        
        %% getR - get model Rs
        function varargout = getR(this, varargin)
            [varargout{1:nargout}] = organism_interface_mex('getR', this.objectHandle, varargin{:});
        end

        %% getParameters - get all model parameters
        function [Parameters] = getParameters(this)
            n = getNParams(this);
            out = cell(n,3);
            for idx = 1:n
              out{idx,1} = idx;
            end
            for idx = 1:n
              out{idx,2} = getParameterName(this,idx);
            end
            for idx = 1:n
              out{idx,3} = getParameterValue(this,idx);
            end
            [Parameters] = out;
        end
    end
end

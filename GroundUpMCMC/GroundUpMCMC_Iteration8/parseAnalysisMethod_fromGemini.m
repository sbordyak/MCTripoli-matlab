function analysisStruct = parseAnalysisMethod(filename)
% parseAnalysisMethod Parses the provided XML file into a MATLAB struct.
%
%   analysisStruct = parseAnalysisMethod(filename) reads the XML file 
%   specified by 'filename' and converts its hierarchical structure into a
%   MATLAB struct. Tags that appear multiple times at the same level (like
%   'ONPEAK') are converted into a struct array.
%
%   Args:
%       filename (string): The path to the XML file.
%
%   Returns:
%       struct: A struct representing the data from the XML file.

    % Check if the file exists
    if ~isfile(filename)
        error('File not found: %s', filename);
    end

    % Read the XML file into a Document Object Model (DOM) node
    try
        domNode = xmlread(filename);
    catch ME
        error('Failed to read or parse XML file: %s\n%s', filename, ME.message);
    end

    % The root element is the first child of the DOM node
    rootElement = domNode.getDocumentElement;

    % Start the recursive parsing from the root element
    analysisStruct = parseChildNodes(rootElement);
end

function dataStruct = parseChildNodes(parentNode)
% parseChildNodes Recursively parses the child nodes of a given DOM node.
%
%   This is a helper function that walks the DOM tree. It converts element
%   nodes into struct fields. If an element contains only text, the text is

%   stored. If it contains other elements, the function calls itself 
%   recursively.
%
%   Args:
%       parentNode (DOM node): The parent node whose children will be parsed.
%
%   Returns:
%       struct: A struct containing the parsed data from the child nodes.
    
    dataStruct = struct();
    childNodes = parentNode.getChildNodes;
    
    for i = 1:childNodes.getLength
        theNode = childNodes.item(i-1);
        
        % Process only element nodes (type 1)
        if theNode.getNodeType == 1
            nodeName = char(theNode.getNodeName);
            
            % Recursively parse the children of the current node
            nodeData = parseChildNodes(theNode);
            
            % If the nodeData is empty, it means this node only has text content
            if isempty(fieldnames(nodeData))
                % Get text content and try to convert to a number
                textContent = strtrim(char(theNode.getTextContent));
                [num, status] = str2num(textContent); %#ok<ST2NM>
                if status
                    nodeData = num;
                else
                    % Handle logicals specifically
                    if strcmpi(textContent, 'true')
                        nodeData = true;
                    elseif strcmpi(textContent, 'false')
                        nodeData = false;
                    else
                        nodeData = textContent;
                    end
                end
            end
            
            % If a field with this name already exists, we need to create a
            % struct array. This handles the multiple <ONPEAK> elements.
            if isfield(dataStruct, nodeName)
                % If it's the first time we see a duplicate, convert the
                % existing single struct into the first element of an array.
                if ~iscell(dataStruct.(nodeName)) && isstruct(dataStruct.(nodeName))
                    dataStruct.(nodeName) = {dataStruct.(nodeName)};
                end
                
                % Append the new struct to the cell array
                dataStruct.(nodeName){end+1} = nodeData;
            else
                % Otherwise, just add the new field
                dataStruct.(nodeName) = nodeData;
            end
        end
    end
    
    % After parsing all children, check if any cell arrays were created.
    % If so, convert them from cell arrays to struct arrays for easier access.
    fields = fieldnames(dataStruct);
    for k = 1:length(fields)
        fieldName = fields{k};
        if iscell(dataStruct.(fieldName))
            try
                % The cat(1, ...) trick concatenates the structs in the
                % cell array into a single struct array.
                dataStruct.(fieldName) = cat(1, dataStruct.(fieldName){:});
            catch ME
                % This might fail if the structs inside the cell are not
                % uniform, in which case we leave it as a cell array.
                warning('Could not convert cell array ''%s'' to a struct array. Leaving as cell. Error: %s', fieldName, ME.message);
            end
        end
    end
    
    % If the top-level struct has only one field (e.g., 'ANALYSIS_METHOD'),
    % it's cleaner to return the contents of that field directly.
    if numel(fields) == 1
        dataStruct = dataStruct.(fields{1});
    end
end

% --- Example Usage ---
% To run this code:
% 1. Save the function above as 'parseAnalysisMethod.m' in your MATLAB path.
% 2. Make sure your XML file 'Pb NBS981 204-5-6-7-8 Daly 14-1-5-5-2 sec.xml' is also in the path.
% 3. Run the following command in the MATLAB command window:
%
%    myAnalysisData = parseAnalysisMethod('Pb NBS981 204-5-6-7-8 Daly 14-1-5-5-2 sec.xml');
%
% 4. You can then inspect the data, for example:
%    disp(myAnalysisData.SETTINGS.TotalBlocks);
%    disp(myAnalysisData.ONPEAK(1).MassID);
%    disp(myAnalysisData.ONPEAK(3).IntegTime);

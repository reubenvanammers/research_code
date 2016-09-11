function pdfsave(handle,name,varargin)
%Saves a image to a pdf, as well as a figure. Places it in the picture
%folder, with resized pages to fit actual figure. Can set options using
%varargin


set(handle,'Units','centimeters','PaperUnits','centimeters',varargin{:});
pos = get(handle,'pos');
set(handle,'PaperPosition',pos,'PaperSize',[pos(3) pos(4)]);
saveas(handle,[pwd '\pictures\' name],'pdf');
savefig(handle,[pwd '\pictures\' name]);


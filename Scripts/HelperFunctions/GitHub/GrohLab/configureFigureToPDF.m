function fig = configureFigureToPDF(fig,figName)
setOpts = {fig,'RendererMode','manual','Renderer','painters',...
        'PaperOrientation','landscape','Color',[1,1,1]};
if exist('figName','var')
    setOpts = cat(2,setOpts,{'Name',figName});
end
set(setOpts{:});
end
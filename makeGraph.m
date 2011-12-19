function makeGraph(name,caption,destdir,relImgDir,xlab,ylab,ylabrule,width,height);
	%This is such a hack, I don't even want to comment on it.
	caption = strrep(caption, '\', '\\');
	wrapper = [
		'\\begin{figure}[htb]\n',...
		'       \\begin{center}\n',...
		'               \\scalebox{0.9}{\n',...
		'                        \\nonstopmode\n',...
		'                        \\input{',relImgDir,'/',name,'.dat.tex}\n',...
		'                        \\errorstopmode\n',...
		'                        \\rule[-0.5cm]{0cm}{0cm}}\n',...
		'                \\caption{',caption,'}\n',...
		'        \\end{center}\n',...
		'\\end{figure}\n'];
	wrapper = strrep(wrapper, '\', '\\');
	wrapper = strrep(wrapper, '$', '\$');

	xlabel(xlab);
	ylabel(['\rule{0pt}{',ylabrule,'}',ylab]);

	system(['mkdir -p ',destdir,' &>/dev/null']);
	print([destdir,'/',name,'.tex'],'-depslatex',['-S',width,',',height]);

	system(['cd ',destdir,'; epstopdf ',name,'.eps; rm ',name,'.eps',...
		"; sed -i 's#",destdir,'#',relImgDir,"#' ", name,".tex"]);

	system(['cd ',destdir,'; ',...
		'mv ',name,'{.tex,.dat.tex}; ',...
		'echo -e "',wrapper,'" > ',name,'.tex;']);
endfunction


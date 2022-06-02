function fig2pdf(filePath)
    % turn all fig files in the input path into pdf with same name, and remove the white border using the pdfcrop tool
    % filePath: absolute path to the folder
    dirOutput = dir(fullfile(filePath, '*.fig'));

    for i = 1:length(dirOutput)
        file_name = sprintf('%s\\%s', dirOutput(i).folder, dirOutput(i).name);
        pdf_file_name = sprintf('%s\\%s', dirOutput(i).folder, replace(dirOutput(i).name, 'fig', 'pdf'));
        fig1 = openfig(file_name);
        print(fig1, '-dpdf', pdf_file_name);
        close;
        status = system(sprintf('pdfcrop %s %s', pdf_file_name, pdf_file_name));
    end

end

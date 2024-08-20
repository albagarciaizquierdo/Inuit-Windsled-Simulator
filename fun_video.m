function fun_video(p,x,filename,t)
    % Crear el objeto de video
    v = VideoWriter(filename, 'MPEG-4');
    v.FrameRate = length(t)/10; % Puedes ajustar la tasa de frames por segundo
    open(v);

    % Iterar sobre cada timestep y dibujar el sistema
    n = size(x, 1);
    box = 0;
    for i = 1:n
        % Crear el título para el frame actual
        title = sprintf('t = %0.1f s', t(i));
        
        % Dibujar el sistema
        fun_draw_system2(p, x(i, :), title, box);
        xlim([-26, 26]);
        ylim([-100, 4]);
        zlim([- 5, 55]);
        
        % Capturar el frame
        frame = getframe(gcf);
        
        % Escribir el frame al video
        writeVideo(v, frame);
        
        % Cerrar la figura
        close(gcf);
    end

    % Cerrar el objeto de video
    close(v);
end

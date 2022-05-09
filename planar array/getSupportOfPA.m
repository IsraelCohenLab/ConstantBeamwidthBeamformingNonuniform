function [S] = getSupportOfPA(M_B, M_A, design)
% S - M_B x M_A boolean support matrix - indicator if there is a sensor at each grid position.
% design options: full, plus, x, star, outline.

S = zeros(M_B, M_A);

if strcmp(design, 'full')
    % All off the sensors are used (M_B * M_A sensors)
    S = ones(M_B, M_A);
    
elseif strcmp(design, 'plus')
    % Only the sensors along the axis are used (M_A + M_B - 1 sensors)
    if mod(M_A,2)==0 || mod(M_B,2)==0
        error('M must be odd');
    end
    
    for n = 1:M_B
        for m = 1:M_A
            if n==(M_B+1)/2 || m==(M_A+1)/2
                S(n,m) = 1;
            end
        end
    end
    
elseif strcmp(design, 'x')
    % Only the sensors along the diagonals are used (M_A + M_B - 1 sensors)
    if M_A~=M_B || mod(M_A,2)==0
        error('M must be odd');
    end
    M = M_A;
    
    for n = 1:M
        for m = 1:M
            if m==n || m==M+1-n
                S(n,m) = 1;
            end
        end
    end
    
elseif strcmp(design, 'star')
    % Generates S for a Star Geometry (based on a linear arrays with M sensors). (4*M - 3 sensors)
    if M_A~=M_B || mod(M_A,2)==0
        error('M must be odd');
    end
    M = M_A;
    
    for n = 1:M
        for m = 1:M
            if m==n || m==M+1-n || n==(M+1)/2 || m==(M+1)/2
                S(n,m) = 1;
            end
        end
    end
    
elseif strcmp(design, 'rectangle_star')
    % Generates S for a Star Geometry (based on a linear arrays with M sensors). (4*M - 3 sensors)
    if mod(M_A,2)==0 || mod(M_B,2)==0
        error('M must be odd');
    end
    N_A = (M_A-1)/2;
    N_B = (M_B-1)/2;
    M = max(M_A, M_B);
    N = (M-1)/2;
    
    for n = -N:N
        for m = -N:N
            if m==n || m==-n || n==0 || m==0
                n_ind = max(min(n,N_B),-N_B)+N_B+1;
                m_ind = max(min(m,N_A),-N_A)+N_A+1;
                S(n_ind,m_ind) = 1;
            end
        end
    end
    
elseif strcmp(design, 'outline')
    % Only the sensors on the perimeter are used (2*M_A + 2*M_B - 4 sensors)
    for n = 1:M_B
        for m = 1:M_A
            if n==1 || n==M_B || m==1 || m==M_A
                S(n,m) = 1;
            end
        end
    end
    
elseif strcmp(design, 'checkered')
    % Sensors are chosen to look like a checkers board.
    % Warning: A_support doesn't seem to have full rank.
    S(1:2:M_B*M_A) = 1; 
    
end
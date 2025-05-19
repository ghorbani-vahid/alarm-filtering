function [attacked_sensor, attack_scenario, ghost_num] = configure_attack(method)
    attack_scenario = input('Ghost attack scenario (deception, replay, delay, none): ', 's');
     if strcmpi(attack_scenario, 'deception')
        ghost_num = input('Number of ghosts: ');
    else
        ghost_num = 0;  % or [] if you prefer an empty indication
    end
    if strcmp(method, 'ST')
        attacked_sensor = 1;
    else
        attacked_sensor = [5, 8];
    end
end
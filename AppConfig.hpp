//
// Created by prana on 25-05-2026.
//

#ifndef APPCONFIG_HPP
#define APPCONFIG_HPP
#include <iostream>
#include <thread>

#endif //APPCONFIG_HPP

#include <fstream>
#include <sstream>

struct AppConfig
{
    int threads = std::thread::hardware_concurrency();
    float noise_floor = 0.08f;
    float attenuation = 0.01f;
};

AppConfig load_config(const std::string& path)
{
    AppConfig cfg;

    std::ifstream file(path);

    if (!file.is_open())
    {
        std::cerr << "Failed to open file " << path << "\n";
        std::cerr << "Proceeding with defaults\n";
        return cfg;

    }

    std::string line;
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }

        std::istringstream ss(line);
        std::string key, value;

        if (!std::getline(ss, key, '=')) continue;
        if (!std::getline(ss, value)) continue;

        key.erase(key.find_last_not_of(" \t") + 1);
        value.erase(0, value.find_first_not_of(" \t"));

        if (key=="threads") cfg.threads = std::stoi(value);
        if (key=="floor") cfg.noise_floor = std::stof(value);
        if (key=="attenuation") cfg.attenuation = std::stof(value);

    }

    std::cerr << "Config loaded from " << path << "\n";
    return cfg;
}
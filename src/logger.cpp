#include "Logger.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>

void Logger::startTimer() {
    t_start = Clock::now();
}

void Logger::stopTimer() {
    t_end = Clock::now();
    durationMs = std::chrono::duration<double, std::milli>(t_end - t_start).count();
}

void Logger::setInputPolygonCount(int count) {
    inputPolygons = count;
}

void Logger::setFaceCount(int count) {
    triangles = count;
}

void Logger::setResultPolygonCount(int count) {
    resultPolygons = count;
}


void Logger::saveToFile(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Fehler beim Öffnen der Log-Datei: " << filename << "\n";
        return;
    }

    out << "===== Ausführungsbericht =====\n";
    out << "Gesamtdauer: " << std::fixed << std::setprecision(3) << durationMs << " ms\n";
    out << "Eingelesene Polygone: " << inputPolygons << "\n";
    out << "Erzeugte Flächen: " << triangles << "\n";
    out << "Ergebnis-Polygone: " << resultPolygons << "\n";

    out.close();
}

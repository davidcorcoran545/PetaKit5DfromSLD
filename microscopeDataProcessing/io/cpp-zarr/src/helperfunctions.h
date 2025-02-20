#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H
#include <string>

#ifndef _WIN32
const char* expandTilde(const char* path);
#endif
    
#ifdef _WIN32
char* strndup (const char *s, size_t n);

int _vscprintf_so(const char * format, va_list pargs);

int vasprintf(char **strp, const char *fmt, va_list ap);

int asprintf(char *strp[], const char *fmt, ...);
#endif

std::string generateUUID();

void mkdirRecursive(const char *dir);

bool fileExists(const std::string &fileName);

void makeDimensionFolders(const std::string &fileName);

bool isLittleEndian();

bool oppositeEndianness(const std::string &dtype);

void swapArrayEndianness(void* array, const size_t elementSize, const size_t numElements);
#endif

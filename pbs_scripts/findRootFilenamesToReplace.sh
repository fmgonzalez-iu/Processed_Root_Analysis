#!/bin/bash
find . -type f | xargs -I{} grep -H -n 'FreeAgentDrive' {}

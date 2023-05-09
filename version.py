#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from os import path, makedirs, unlink
from datetime import datetime, date
from zipfile import ZipFile

parser = argparse.ArgumentParser()
parser.add_argument('todo',
                    nargs    = 1,
                    type     = str,
                    metavar  = 'TODO',
                    help     = "'standalone' makes standalone zip; major, minor, patch bump version.")
parser.add_argument('--dry',
                    dest     = 'dry',
                    action   = 'store_true',
                    help     = "Dry run (do not run)",
                    required = False)
args = vars(parser.parse_args())

# ---------------------------------------------------------------------
# Config

config_version = "0.6.8"
config_date = date(2023, 5, 9)
config_files = [
    ('version.py', 'config_version = "{major}.{minor}.{patch}"'),
    ('version.py', 'config_date = date({date:%Y, %-m, %-d})'),
    ('README.md', 'version {major}.{minor}.{patch} {date:%d%b%Y}'),
    ('manyiv.pkg', 'v {major}.{minor}.{patch}'),
    ('manyiv.pkg', 'd Distribution-Date: {date:%Y%m%d}'),
    ('stata.toc', 'v {major}.{minor}.{patch}'),
    ('src/ado/manyiv.sthlp', 'version {major}.{minor}.{patch} {date:%d%b%Y}'),
    ('src/ado/manyiv.ado', 'version {major}.{minor}.{patch} {date:%d%b%Y}')
]

config_standalone = {
    'manyiv': [
        'src/build/lmanyiv.mlib',
        'src/ado/manyiv.ado',
        'src/ado/manyiv.sthlp',
        'src/mata/manyiv_internals.mata',
        'src/mata/manyiv_absorb.mata',
        'src/build/manyiv_unix.plugin',
        'src/build/manyiv_macosx.plugin',
        'src/build/manyiv_windows.plugin'
    ]
}

# ---------------------------------------------------------------------
# Bump or standalone


def main(todo, dry = False):
    args = ['major', 'minor', 'patch', 'standalone']
    if todo not in args:
        msg = f"'{todo}' uknown; can only bump: {', '.join(args)}"
        raise Warning(msg)

    if todo == 'standalone':
        make_standalone(config_standalone, dry)
    else:
        current_kwargs, update_kwargs = bump_kwargs(todo, config_version, config_date)
        for file, string in config_files:
            bump_file(file, string, current_kwargs, update_kwargs, dry)


def make_standalone(standalone, dry = False):
    for label, files in standalone.items():
        outzip = f'standalone/{label}-{config_version}.zip'

        if dry:
            for f in files:
                print(f'{f} -> {outzip}')
        else:
            if not path.isdir('standalone'):
                makedirs('standalone')

            if path.isfile(outzip):
                unlink(outzip)

            with ZipFile(outzip, 'w') as zf:
                for f in files:
                    print(f'{f} -> {outzip}')
                    zf.write(f, path.basename(f))


def bump_file(file, string, current, update, dry = False):
    find = string.format(**current)
    with open(file, 'r') as fh:
        lines = fh.readlines()
        if find not in ''.join(lines):
            print(f'WARNING: nothing to bump in {file}')

        replace = string.format(**update)
        ulines = []
        for line in lines:
            if find in line:
                print(f'{file}: {find} -> {replace}')
                ulines += [line.replace(find, replace)]
            else:
                ulines += [line]

    if not dry:
        with open(file, 'w') as fh:
            fh.write(''.join(ulines))


def bump_kwargs(bump, config_version, config_date):
    today = datetime.now()
    major, minor, patch = config_version.split('.')
    umajor, uminor, upatch = bump_sever(bump, major, minor, patch)

    current_kwargs = {
        'major': major,
        'minor': minor,
        'patch': patch,
        'date': config_date
    }

    update_kwargs = {
        'major': umajor,
        'minor': uminor,
        'patch': upatch,
        'date': today
    }

    return current_kwargs, update_kwargs


def bump_sever(bump, major, minor, patch):
    if bump == 'major':
        return str(int(major) + 1), '0', '0'
    elif bump == 'minor':
        return major, str(int(minor) + 1), '0'
    elif bump == 'patch':
        return major, minor, str(int(patch) + 1)
    else:
        return major, minor, patch


if __name__ == "__main__":
    main(args['todo'][0], args['dry'])

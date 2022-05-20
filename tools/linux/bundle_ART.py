#!/usr/bin/env python3

import os, sys
import shutil
import subprocess
import argparse

def getdlls(opts):
    blacklist = {
	'linux-vdso.so.1',
	'libm.so.6',
	'libpthread.so.0',
	'libc.so.6',
	'ld-linux-x86-64.so.2',
	'libdl.so.2',
	'libsystemd.so.0',
        'librt.so.1',
        'libstdc++.so.6',
        }
    res = []
    d = os.getcwd()
    p = subprocess.Popen(['ldd', os.path.join(d, 'ART')],
                         stdout=subprocess.PIPE)
    out, _ = p.communicate()
    for line in out.decode('utf-8').splitlines():
        if ' => ' in line:
            bits = line.split()
            lib = bits[2]
            bn = os.path.basename(lib)
            if bn not in blacklist:
                res.append(lib)
    return res

def getopts():
    p = argparse.ArgumentParser()
    p.add_argument('-o', '--outdir', required=True,
                   help='output directory for the bundle')
    p.add_argument('-e', '--exiftool', help='path to exiftool dir')
    p.add_argument('-v', '--verbose', action='store_true')
    ret = p.parse_args()
    return ret

def extra_files(opts):
    def D(s): return os.path.expanduser(s)
    if opts.exiftool and os.path.isdir(opts.exiftool):
        extra = [('lib', [(opts.exiftool, 'exiftool')])]
    else:
        extra = []
    return [
        ('share/icons/Adwaita', [
            D('/usr/share/icons/Adwaita/scalable'),
            D('/usr/share/icons/Adwaita/index.theme'), 
            D('/usr/share/icons/Adwaita/cursors'),
        ]),
        ('lib', [
            D('/usr/lib/x86_64-linux-gnu/gdk-pixbuf-2.0'),
        ]),
        ('lib', [
            D('/usr/lib/x86_64-linux-gnu/gio'),
        ]),
        ('share/glib-2.0/schemas', [
            D('/usr/share/glib-2.0/schemas/gschemas.compiled'),
        ]),
        ('share', [
            (D('~/.local/share/lensfun/updates/version_1'), 'lensfun'),
        ]),
        ('lib', [
            D('/usr/lib/x86_64-linux-gnu/gvfs/libgvfscommon.so'),
            D('/usr/lib/x86_64-linux-gnu/gvfs/libgvfsdaemon.so'),
        ]),
    ] + extra


def main():
    opts = getopts()
    d = os.getcwd()
    if not os.path.exists('ART'):
        sys.stderr.write('ERROR: ART not found! Please run this script '
                         'from the build directory of ART\n')
        sys.exit(1)
    if opts.verbose:
        print('copying %s to %s' % (os.getcwd(), opts.outdir))
    shutil.copytree(d, opts.outdir)
    if not os.path.exists(os.path.join(opts.outdir, 'lib')):
        os.mkdir(os.path.join(opts.outdir, 'lib'))
    for lib in getdlls(opts):
        if opts.verbose:
            print('copying: %s' % lib)
        shutil.copy2(lib,
                     os.path.join(opts.outdir, 'lib', os.path.basename(lib)))
    for key, elems in extra_files(opts):
        for elem in elems:
            name = None
            if isinstance(elem, tuple):
                elem, name = elem
            else:
                name = os.path.basename(elem)
            if opts.verbose:
                print('copying: %s' % elem)
            if not os.path.exists(elem):
                print('SKIPPING non-existing: %s' % elem)
            elif os.path.isdir(elem):
                shutil.copytree(elem, os.path.join(opts.outdir, key, name))
            else:
                dest = os.path.join(opts.outdir, key, name)
                destdir = os.path.dirname(dest)
                if not os.path.exists(destdir):
                    os.makedirs(destdir)
                shutil.copy2(elem, dest)
    os.makedirs(os.path.join(opts.outdir, 'share/gtk-3.0'))
    with open(os.path.join(opts.outdir, 'share/gtk-3.0/settings.ini'), 'w') \
         as out:
        out.write('[Settings]\ngtk-button-images=1\n')
    with open(os.path.join(opts.outdir, 'options'), 'a') as out:
        out.write('\n[Lensfun]\nDBDirectory=share/lensfun\n')
        if opts.exiftool:
            out.write('\n[Metadata]\nExiftoolPath=exiftool\n')
    for name in ('ART', 'ART-cli'):
        shutil.move(os.path.join(opts.outdir, name),
                    os.path.join(opts.outdir, name + '.bin'))
    with open(os.path.join(opts.outdir, 'fonts.conf'), 'w') as out:
        out.write("""\
<?xml version="1.0"?>
<!DOCTYPE fontconfig SYSTEM "fonts.dtd">

<fontconfig>
  <dir>/usr/share/fonts</dir>
  <dir>/usr/local/share/fonts</dir>
  <dir prefix="xdg">fonts</dir>

  <cachedir>/var/cache/fontconfig</cachedir>
  <cachedir prefix="xdg">fontconfig</cachedir>

  <match target="pattern">
    <test qual="any" name="family"><string>mono</string></test>
    <edit name="family" mode="assign" binding="same"><string>monospace</string></edit>
  </match>
  <match target="pattern">
    <test qual="any" name="family"><string>sans serif</string></test>
    <edit name="family" mode="assign" binding="same"><string>sans-serif</string></edit>
  </match>

  <match target="pattern">
    <test qual="any" name="family"><string>sans</string></test>
    <edit name="family" mode="assign" binding="same"><string>sans-serif</string></edit>
  </match>

  <alias>
    <family>DejaVu Sans</family>
    <default><family>sans-serif</family></default>
  </alias>
  <alias>
    <family>sans-serif</family>
    <prefer><family>DejaVu Sans</family></prefer>
  </alias>
  <alias>
    <family>DejaVu Serif</family>
    <default><family>serif</family></default>
  </alias>
  <alias>
    <family>serif</family>
    <prefer><family>DejaVu Serif</family></prefer>
  </alias>
  <alias>
    <family>DejaVu Sans Mono</family>
    <default><family>monospace</family></default>
  </alias>
  <alias>
    <family>monospace</family>
    <prefer><family>DejaVu Sans Mono</family></prefer>
  </alias>

  <config><rescan><int>30</int></rescan></config>
</fontconfig>
""")        
    with open(os.path.join(opts.outdir, 'ART'), 'w') as out:
        out.write("""#!/bin/bash
export ART_restore_GTK_CSD=$GTK_CSD
export ART_restore_GDK_PIXBUF_MODULE_FILE=$GDK_PIXBUF_MODULE_FILE
export ART_restore_GDK_PIXBUF_MODULEDIR=$GDK_PIXBUF_MODULEDIR
export ART_restore_GIO_MODULE_DIR=$GIO_MODULE_DIR
export ART_restore_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
export ART_restore_FONTCONFIG_FILE=$FONTCONFIG_FILE
export GTK_CSD=0
d=$(dirname $(readlink -f "$0"))
t=$(mktemp -d --suffix=-ART)
"$d/lib/gdk-pixbuf-2.0/gdk-pixbuf-query-loaders" "$d/lib/gdk-pixbuf-2.0/2.10.0/loaders/libpixbufloader-png.so" "$d/lib/gdk-pixbuf-2.0/2.10.0/loaders/libpixbufloader-svg.so" > "$t/loader.cache"
export GDK_PIXBUF_MODULE_FILE="$t/loader.cache"
export GDK_PIXBUF_MODULEDIR="$d/lib/gdk-pixbuf-2.0"
export GIO_MODULE_DIR="$d/lib/gio/modules"
export LD_LIBRARY_PATH="$d/lib"
export FONTCONFIG_FILE="$d/fonts.conf"
export ART_EXIFTOOL_BASE_DIR="$d/lib/exiftool"
"$d/ART.bin" "$@"
rm -rf "$t"
""")
    with open(os.path.join(opts.outdir, 'ART-cli'), 'w') as out:
        out.write("""#!/bin/bash
export ART_restore_GIO_MODULE_DIR=$GIO_MODULE_DIR
export ART_restore_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
d=$(dirname $(readlink -f "$0"))
export GIO_MODULE_DIR="$d/lib/gio/modules"
export LD_LIBRARY_PATH="$d/lib"
export ART_EXIFTOOL_BASE_DIR="$d/lib/exiftool"
exec "$d/ART-cli.bin" "$@"
""")
    for name in ('ART', 'ART-cli'):
        os.chmod(os.path.join(opts.outdir, name), 0o755)

if __name__ == '__main__':
    main()

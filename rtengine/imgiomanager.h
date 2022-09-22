/* -*- C++ -*-
 *
 *  This file is part of ART.
 *
 *  Copyright 2020 Alberto Griggio <alberto.griggio@gmail.com>
 *
 *  ART is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ART is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ART.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "noncopyable.h"
#include "rtengine.h"
#include "imageio.h"
#include <glibmm/ustring.h>
#include <unordered_map>
#include <map>

namespace rtengine {

class ImageIOManager: public NonCopyable {
public:
    enum Format {
        FMT_UNKNOWN = 0,
        FMT_JPG,
        FMT_PNG,
        FMT_PNG16,
        FMT_TIFF,
        FMT_TIFF_FLOAT,
        FMT_TIFF_FLOAT16
    };

    struct SaveFormatInfo {
        std::string extension;
        Glib::ustring label;

        SaveFormatInfo(const std::string &ext="", const Glib::ustring &lbl=""):
            extension(ext), label(lbl) {}
    };

    static ImageIOManager *getInstance();

    void init(const Glib::ustring &dir);
    bool load(const Glib::ustring &fileName, ProgressListener *plistener, ImageIO *&img, int maxw_hint, int maxh_hint);
    bool save(IImagefloat *img, const std::string &ext, const Glib::ustring &fileName, ProgressListener *plistener);
    std::vector<std::pair<std::string, SaveFormatInfo>> getSaveFormats() const;

    bool canLoad(const std::string &ext) const
    {
        return loaders_.find(ext) != loaders_.end();
    }

    Format getFormat(const Glib::ustring &fileName);

private:
    static Glib::ustring get_ext(Format f);
    
    std::unordered_map<std::string, Glib::ustring> loaders_;
    std::unordered_map<std::string, Glib::ustring> savers_;
    std::unordered_map<std::string, Format> fmts_;
    std::map<std::string, SaveFormatInfo> savelbls_;
    
    Glib::ustring dir_;
};

} // namespace rtengine

// [[file:../edip.note::1141962a][1141962a]]
//#![deny(warnings)]
mod edip;
// 1141962a ends here

// [[file:../edip.note::*mods][mods:1]]

// mods:1 ends here

// [[file:../edip.note::3d96dcdb][3d96dcdb]]
pub use crate::edip::*;
// 3d96dcdb ends here

// [[file:../edip.note::*docs][docs:1]]
#[cfg(feature = "adhoc")]
/// Docs for local mods
pub mod docs {
    macro_rules! export_doc {
        ($l:ident) => {
            pub mod $l {
                pub use crate::$l::*;
            }
        };
    }

    // export_doc!(codec);
}
// docs:1 ends here

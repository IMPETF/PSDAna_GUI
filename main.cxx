#include "TApplication.h"
#include "GuiFrame.h"
#include "TGClient.h"

int main(int argc, char *argv[])
{

    TApplication apponline("Online", &argc, argv);

    GuiFrame *mygui = new GuiFrame(gClient->GetRoot(), 250, 250);

    apponline.Run();

    return 0;
}
